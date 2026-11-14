/*
 * calc.js -- the edms calculators, in the browser.
 *
 * These mirror the Python in edms.work.plates and edms.bio.*, so a rendered
 * protocol can recalculate without a server. That duplication is the whole
 * risk of this file: a number that drifts from the Python is worse than no
 * number at all, because it looks authoritative. Two things hold it in place.
 *
 * 1. Every constant comes from Python at render time, in EDMS.consts. Nothing
 *    numeric is written down twice.
 * 2. tests/work/test_calc_js.py runs this file under node and compares its
 *    output to the Python's, over the protocols' own values and over
 *    randomised ones. A port that drifts fails the suite.
 *
 * Anything not ported here is not guessed at: the engine marks the table
 * static and the CLI stays the way to recompute it.
 */
'use strict';

var EDMS = (typeof EDMS !== 'undefined') ? EDMS : {};
EDMS.calc = {};

(function (C) {

  // ----------------------------------------------------------------- //
  // Numbers
  // ----------------------------------------------------------------- //
  /* Python's round() is half-to-even where JavaScript's Math.round() is
     half-up, so round(0.125, 2) is 0.12 there and would be 0.13 here. That
     much is well known. The part that actually bites is subtler: rounding by
     multiplying by a power of ten destroys the information the decision needs.
     4.4 * 1.1875 is 5.2250000000000005, which Python rounds up to 5.23, but
     times 100 it lands on exactly 522.5 and a half-to-even rule then rounds it
     down to 5.22. The reverse happens to 2.675, which is really 2.67499... and
     must round down, yet times 100 becomes exactly 267.5 and rounds up.

     So work from the double's own decimal expansion, which toFixed gives
     exactly, and decide from the digits that are actually there. Every
     rounding in a ported calculator goes through this. */
  C.pyRound = function (x, n) {
    n = n || 0;
    if (x === null || x === undefined || typeof x !== 'number' || !isFinite(x))
      return x;
    if (x === 0) return 0;
    // Outside this range the digit string stops fitting a safe integer; no
    // lab volume gets there, and the plain form is right for all of them.
    if (n < 0 || n > 15 || !(Math.abs(x) * Math.pow(10, n) < 1e15)) {
      var f0 = Math.pow(10, n);
      return Math.round(x * f0) / f0;
    }
    var neg = x < 0;
    var s = Math.abs(x).toFixed(Math.min(100, n + 25));
    var dot = s.indexOf('.');
    var head = s.slice(0, dot) + s.slice(dot + 1, dot + 1 + n);
    var rest = s.slice(dot + 1 + n);
    var up;
    var first = rest.charAt(0);
    if (first > '5') up = true;
    else if (first < '5') up = false;
    else if (/[1-9]/.test(rest.slice(1))) up = true;  // just over the half
    else up = (parseInt(head.charAt(head.length - 1), 10) % 2) === 1;  // to even
    var val = parseInt(head, 10) + (up ? 1 : 0);
    var out = val / Math.pow(10, n);
    return neg ? -out : out;
  };

  /* Python's round() returns a float, so it prints 10.0 where String(10)
     here would print 10. These go into prose cells that are compared as
     text, so the trailing .0 is not cosmetic. */
  C.pyFloat = function (v) {
    return Number.isInteger(v) ? v.toFixed(1) : String(v);
  };

  C.num = function (v) {
    if (v === null || v === undefined || v === '' || v === 'N/A' || v === '-')
      return null;
    var x = typeof v === 'number' ? v : parseFloat(v);
    return isFinite(x) ? x : null;
  };

  // ----------------------------------------------------------------- //
  // Formats (edms.work.plates)
  // ----------------------------------------------------------------- //
  C.canonical = function (fmt) {
    var key = String(fmt).trim().toLowerCase().split(/\s+/).join(' ');
    var a = EDMS.consts.aliases;
    if (Object.prototype.hasOwnProperty.call(a, key)) return a[key];
    throw new Error("Unknown plate format '" + fmt + "'");
  };

  C.spec = function (fmt) {
    var name = C.canonical(fmt);
    var v = EDMS.consts.vessels[name];
    return {plate_format: name, wells: v.wells, area_cm2: v.area_cm2,
            media_mL: v.media_mL,
            seed_HEK293T: C.roundCells(v.area_cm2 * EDMS.consts.SEED_PER_CM2)};
  };

  C.scale = function (to, frm) {
    return C.spec(to).area_cm2 / C.spec(frm).area_cm2;
  };

  C.names = function () {
    var ns = Object.keys(EDMS.consts.vessels);
    ns.sort(function (a, b) {
      return EDMS.consts.vessels[a].area_cm2 - EDMS.consts.vessels[b].area_cm2;
    });
    return ns;
  };

  /* mirrors plates._round_cells */
  C.roundCells = function (cells) {
    if (cells >= 1e5) {
      var step = Math.pow(10, Math.floor(Math.log10(cells)) - 1);
      return Math.round(C.pyRound(cells / step, 0) * step);
    }
    if (cells >= 1e4) return Math.round(C.pyRound(cells / 1000, 0) * 1000);
    return Math.round(C.pyRound(cells / 100, 0) * 100);
  };

  // ----------------------------------------------------------------- //
  // plates:geometry
  // ----------------------------------------------------------------- //
  C['edms.work.plates:geometry'] = function (a) {
    var s = C.spec(a.plate_format);
    var reps = parseInt(a.reps === undefined ? 1 : a.reps, 10);
    var density = (a.seed_per_cm2 === undefined || a.seed_per_cm2 === null)
      ? EDMS.consts.SEED_PER_CM2 : parseFloat(a.seed_per_cm2);
    var seed = C.roundCells(s.area_cm2 * density);
    var row = {
      'Format': s.plate_format,
      'Growth area (cm2)': s.area_cm2,
      'Media per well (mL)': s.media_mL,
      'Seed HEK293T per well': seed,
      'Seeding density (per cm2)': Math.round(C.pyRound(density, 0)),
      'Replicates': reps
    };
    var cols = Object.keys(row);
    if (a.conditions !== undefined && a.conditions !== null && a.conditions !== '') {
      var conditions = parseInt(a.conditions, 10);
      var wells = conditions * reps;
      row['Conditions'] = conditions;
      row['Wells needed'] = wells;
      row['Plates needed'] = Math.ceil(wells / s.wells);
      row['Media needed (mL)'] = C.pyRound(wells * s.media_mL, 1);
      row['Cells needed'] = wells * seed;
      cols = Object.keys(row);
    }
    return {columns: cols, rows: [row]};
  };

  // ----------------------------------------------------------------- //
  // plates:virus  (wraps edms.bio.transfect:virus)
  // ----------------------------------------------------------------- //
  C['edms.work.plates:virus'] = function (a) {
    var f = C.scale(a.plate_format || EDMS.consts.VIRUS_ANCHOR,
                    EDMS.consts.VIRUS_ANCHOR);
    var d = function (v, dflt) {
      return C.pyRound((v === undefined || v === null ? dflt : parseFloat(v)) * f, 4);
    };
    return C._transfectVirus({
      plasmids: a.plasmids,
      VSVG_plasmid: a.VSVG_plasmid || 'pMUZ26.6',
      GagPol_plasmid: a.GagPol_plasmid || 'pMUZ26.7',
      reps: parseInt(a.reps === undefined ? 1 : a.reps, 10),
      mm_x: a.mm_x === undefined ? 1.1 : parseFloat(a.mm_x),
      VSVG_ng: d(a.VSVG_ng, 750),
      GagPol_ng: d(a.GagPol_ng, 1500),
      transfer_ng: d(a.transfer_ng, 750),
      well_uL: d(a.well_uL, 500),
      P3000_uL_per_well: d(a.P3000_uL_per_well, 6),
      L3000_uL_per_well: d(a.L3000_uL_per_well, 7)
    });
  };

  C._transfectVirus = function (a) {
    var rows = a.plasmids.rows || a.plasmids;
    var transfers = rows.filter(function (r) {
      return r.Plasmid !== a.VSVG_plasmid && r.Plasmid !== a.GagPol_plasmid;
    }).length;
    var batch = transfers * a.mm_x * a.reps;
    var optimem = C.pyRound(batch * a.well_uL / 2, 2);
    var p3000 = C.pyRound(batch * a.P3000_uL_per_well, 2);
    var l3000 = C.pyRound(batch * a.L3000_uL_per_well, 2);

    var body = rows.map(function (r) {
      var conc = C.num(r['ng/uL']);
      var ng, uL, tube, order;
      if (r.Plasmid === a.VSVG_plasmid) {
        ng = C.pyRound(batch * a.VSVG_ng, 2);
        uL = C.pyRound(batch * a.VSVG_ng / conc, 2);
        tube = 'A'; order = 1;
      } else if (r.Plasmid === a.GagPol_plasmid) {
        ng = C.pyRound(batch * a.GagPol_ng, 2);
        uL = C.pyRound(batch * a.GagPol_ng / conc, 2);
        tube = 'A'; order = 1;
      } else {
        ng = C.pyRound(a.transfer_ng * a.reps, 2);
        uL = C.pyRound(a.transfer_ng * a.reps / conc, 2);
        tube = (String(r.Plasmid).match(/\d+\.\d+|\d+/g) || []).join(', ');
        order = 2;
      }
      return {Plasmid: r.Plasmid, Description: r.Description,
              Colony: r.Colony, 'ng/uL': conc, ng: ng, uL: uL,
              Tube: tube, Order: order};
    });
    // pandas sort_values(by=['Order', 'Plasmid']), which is stable
    body.sort(function (x, y) {
      if (x.Order !== y.Order) return x.Order - y.Order;
      return x.Plasmid < y.Plasmid ? -1 : (x.Plasmid > y.Plasmid ? 1 : 0);
    });

    var pad = function (desc, uL, tube, order) {
      return {Plasmid: 'N/A', Description: desc, Colony: 'N/A', 'ng/uL': 'N/A',
              ng: 'N/A', uL: uL, Tube: tube, Order: order};
    };
    var out = [pad('Optimem', optimem, 'A', 1), pad('P3000 Enhancer', p3000, 'A', 1)]
      .concat(body, [pad('Optimem', optimem, 'B', 3), pad('L3000', l3000, 'B', 3)]);
    return {columns: ['Plasmid', 'Description', 'Colony', 'ng/uL', 'ng', 'uL',
                      'Tube', 'Order'], rows: out};
  };

  // ----------------------------------------------------------------- //
  // plates:virus_per_well
  // ----------------------------------------------------------------- //
  C['edms.work.plates:virus_per_well'] = function (a) {
    var fmt = a.plate_format || EDMS.consts.VIRUS_ANCHOR;
    var f = C.scale(fmt, EDMS.consts.VIRUS_ANCHOR);
    var s = C.spec(fmt);
    var g = function (v, dflt) {
      return (v === undefined || v === null || v === '') ? dflt : parseFloat(v);
    };
    var wellUL = g(a.well_uL, 500);
    var dna = g(a.VSVG_ng, 750) + g(a.GagPol_ng, 1500) + g(a.transfer_ng, 750);
    var complexUL = wellUL * f;
    var row = {
      'Format': s.plate_format,
      'Growth area (cm2)': s.area_cm2,
      'Complex per well (uL)': C.pyRound(complexUL, 2),
      'Media per well (mL)': s.media_mL,
      'Complex as % of media': C.pyRound(100 * complexUL / (s.media_mL * 1000), 1),
      'DNA per well (ng)': C.pyRound(dna * f, 1),
      'DNA per cm2 (ng)': C.pyRound(dna * f / s.area_cm2, 1),
      'P3000 per well (uL)': C.pyRound(g(a.P3000_uL_per_well, 6) * f, 3),
      'L3000 per well (uL)': C.pyRound(g(a.L3000_uL_per_well, 7) * f, 3)
    };
    return {columns: Object.keys(row), rows: [row]};
  };

  // ----------------------------------------------------------------- //
  // plates:dilutions / virus_dilutions
  // ----------------------------------------------------------------- //
  var DIL_COLS = ['Component', 'Planned uL', 'Stock ng/uL', 'Fold',
                  'Working ng/uL', 'uL at working', 'Make'];

  C.dilutions = function (plan, a) {
    var minUL = a.min_uL === undefined || a.min_uL === null || a.min_uL === ''
      ? EDMS.consts.MIN_PIPETTE_UL : parseFloat(a.min_uL);
    var makeUL = a.make_uL === undefined || a.make_uL === null || a.make_uL === ''
      ? 20.0 : parseFloat(a.make_uL);
    var diluent = a.diluent === undefined ? 'TE' : a.diluent;
    if (!(minUL > 0)) throw new Error('min_uL must be greater than 0');

    var out = [];
    (plan.rows || plan).forEach(function (r) {
      var uL = C.num(r.uL);
      if (uL === null || uL >= minUL || uL <= 0) return;
      var name = (r.Plasmid === undefined || String(r.Plasmid) === 'N/A')
        ? r.Description : r.Plasmid;
      var conc = C.num(r['ng/uL']);
      var ng = C.num(r.ng);
      var exact = (conc && ng !== null && conc > 0) ? ng / conc : uL;
      var needed = minUL / exact;
      if (conc === null || conc <= 0 || ng === null) {
        out.push({'Component': name, 'Planned uL': C.pyRound(uL, 3),
                  'Stock ng/uL': '-', 'Fold': '-', 'Working ng/uL': '-',
                  'uL at working': '-',
                  'Make': 'scale the mix ' + needed.toFixed(1) +
                          'x: plan more wells, or raise mm_x'});
        return;
      }
      var fold = null;
      for (var i = 0; i < EDMS.consts.DILUTION_STEPS.length; i++) {
        if (EDMS.consts.DILUTION_STEPS[i] >= needed) {
          fold = EDMS.consts.DILUTION_STEPS[i];
          break;
        }
      }
      if (fold === null) fold = Math.ceil(needed / 1000) * 1000;
      var working = conc / fold;
      var total = Math.max(makeUL, minUL * fold);
      var stockUL = total / fold;
      out.push({'Component': name, 'Planned uL': C.pyRound(uL, 3),
                'Stock ng/uL': C.pyRound(conc, 2), 'Fold': fold,
                'Working ng/uL': C.pyRound(working, 2),
                'uL at working': C.pyRound(ng / working, 2),
                'Make': C.pyFloat(C.pyRound(stockUL, 2)) + ' uL stock + ' +
                        C.pyFloat(C.pyRound(total - stockUL, 2)) + ' uL ' +
                        diluent});
    });
    return {columns: DIL_COLS, rows: out};
  };

  C['edms.work.plates:virus_dilutions'] = function (a) {
    return C.dilutions(C['edms.work.plates:virus'](a), a);
  };

  // ----------------------------------------------------------------- //
  // plates:expansion
  // ----------------------------------------------------------------- //
  C['edms.work.plates:expansion'] = function (a) {
    var order = C.names();
    var first = C.canonical(a.start), last = C.canonical(a.end);
    if (order.indexOf(last) <= order.indexOf(first))
      throw new Error("'" + last + "' is not larger than '" + first + "'");
    var skipped = (a.skip || EDMS.consts.EXPANSION_SKIP_DEFAULT || [])
      .map(function (s) { return C.canonical(s); });
    var rungs = order.slice(order.indexOf(first), order.indexOf(last) + 1)
      .filter(function (n) { return skipped.indexOf(n) === -1; });
    var rows = [];
    for (var i = 0; i + 1 < rungs.length; i++) {
      var frm = rungs[i], to = rungs[i + 1];
      var ratio = C.scale(to, frm);
      var between = order.slice(order.indexOf(frm) + 1, order.indexOf(to))
        .filter(function (n) { return skipped.indexOf(n) !== -1; });
      var note = between.length ? 'skipping ' + between.join(', ') : '';
      if (ratio >= EDMS.consts.WIDE_SPLIT) {
        var wide = '1:' + C.pyRound(ratio, 0).toFixed(0) + ' is a wide split; ' +
          'expect to sit at low density for a few days';
        note = note ? note + ', ' + wide : wide;
      }
      rows.push({'From': frm, 'To': to, 'Split': '1:' + ratio.toFixed(1),
                 'Area ratio': C.pyRound(ratio, 2),
                 'Media per well (mL)': C.spec(to).media_mL, 'Note': note});
    }
    return {columns: ['From', 'To', 'Split', 'Area ratio',
                      'Media per well (mL)', 'Note'], rows: rows};
  };

})(EDMS.calc);

if (typeof module !== 'undefined' && module.exports) module.exports = EDMS;
