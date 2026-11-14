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

  /* Python's format(x, 'g'): six significant digits, trailing zeros stripped,
     and an exponent once the value gets far enough from 1. It shows up in the
     prose cells -- '100 ng at 412 ng/uL' -- where 100.0 would read wrong. */
  C.pyG = function (x) {
    if (typeof x !== 'number') return String(x);
    if (!isFinite(x)) return x > 0 ? 'inf' : (x < 0 ? '-inf' : 'nan');
    if (x === 0) return '0';
    var r = Math.abs(x).toExponential(5);
    // Read the exponent back off the rounded form: 999999.9 rounds to 1e+06
    var exp = parseInt(r.slice(r.indexOf('e') + 1), 10);
    var out;
    if (exp < -4 || exp >= 6) {
      var m = r.slice(0, r.indexOf('e'));
      if (m.indexOf('.') !== -1) m = m.replace(/0+$/, '').replace(/\.$/, '');
      out = m + 'e' + (exp < 0 ? '-' : '+') +
        (Math.abs(exp) < 10 ? '0' : '') + Math.abs(exp);
    } else {
      out = Math.abs(x).toFixed(Math.max(0, 5 - exp));
      if (out.indexOf('.') !== -1) out = out.replace(/0+$/, '').replace(/\.$/, '');
    }
    return (x < 0 ? '-' : '') + out;
  };

  /* A value printed into a label the way Python's f-string would print it.
     JSON has one number and Python has two, so a template that writes
     mm_x = 1.0 wants '1.0' in the master mix heading where String() gives
     '1'. EDMS.floats, threaded through resolveArgs, carries the difference;
     a value the user typed is already a string and goes through as typed,
     because edms work deviate parses it back to exactly that. */
  C.label = function (a, name, dflt) {
    var v = (a && a[name] !== undefined && a[name] !== null && a[name] !== '')
      ? a[name] : dflt;
    if (typeof v === 'string') return v;
    if (typeof v !== 'number') return String(v);
    return ((a.__floats__ || []).indexOf(name) !== -1 || !Number.isInteger(v))
      ? C.pyFloat(v) : String(v);
  };

  /* The value a table echoes back, as opposed to the one it calculates
     with. Python puts the parameter itself in the cell, so whatever came in
     goes out: a number stays a number and a string the user typed stays the
     string, which is what edms work deviate would have parsed it back to. */
  C.raw = function (a, name, dflt) {
    var v = a[name];
    return (v === undefined || v === null || v === '') ? dflt : v;
  };

  // ----------------------------------------------------------------- //
  // Frames
  // ----------------------------------------------------------------- //
  /* A calculator that takes a table gets one of two things: a payload table
     read out of the run's inputs, or rows this file built itself. pandas
     turns a list of dicts into a frame whose columns are the keys in the
     order they first appear, and a key one row omits reads as NaN there. The
     Python asks 'is this column in the frame at all' separately from 'does
     this row have a value', and the two answers differ -- an absent 'In MM'
     column means every component is in the mix. Keep the distinction. */
  C.frame = function (obj) {
    if (obj && obj.columns && obj.rows) return obj;
    var rows = obj || [];
    var columns = [];
    rows.forEach(function (r) {
      Object.keys(r).forEach(function (k) {
        if (columns.indexOf(k) === -1) columns.push(k);
      });
    });
    return {columns: columns, rows: rows};
  };

  C.has = function (f, col) { return f.columns.indexOf(col) !== -1; };

  /* A cell's value, or null where Python would see NaN, an empty string or no
     such column. */
  C.cell = function (f, row, col) {
    if (!C.has(f, col)) return null;
    var v = row[col];
    return (v === null || v === undefined || v === '') ? null : v;
  };

  C.cellNum = function (f, row, col) {
    var v = C.cell(f, row, col);
    if (v === null) return null;
    var x = typeof v === 'number' ? v : parseFloat(v);
    return isFinite(x) ? x : null;
  };

  C.blank = function (v) { return (v === null || v === undefined) ? '' : v; };

  C.g = function (a, name, dflt) {
    var v = a[name];
    return (v === undefined || v === null || v === '') ? dflt : parseFloat(v);
  };

  // ----------------------------------------------------------------- //
  // clone:rxn_mm -- the engine under every reaction table
  // ----------------------------------------------------------------- //
  C._rxnMM = function (input, o) {
    var f = C.frame(input);
    var cCol = o.component_col || 'Component';
    var sCol = o.stock_col || 'Stock';
    var dCol = o.desired_col || 'Desired';
    var uCol = o.unit_col || 'Unit';
    var vCol = o.uL_col || 'uL';
    var mCol = o.mm_col || 'In MM';
    var nCol = o.note_col || 'Note';
    var water = o.water || 'Nuclease-free water';
    var rxn = parseFloat(o.rxn_uL);
    var samples = parseInt(o.samples === undefined ? 1 : o.samples, 10);
    var mmX = (o.mm_x === undefined || o.mm_x === null || o.mm_x === '')
      ? 1.1 : parseFloat(o.mm_x);

    if (!C.has(f, cCol))
      throw new Error("rxn_mm() needs a '" + cCol + "' column; df has " +
                      f.columns.join(', '));
    if (!(rxn > 0))
      throw new Error('rxn_uL must be greater than 0, got ' + rxn);
    if (!(samples >= 1))
      throw new Error('samples must be at least 1, got ' + samples);

    var names = [], stocks = [], desireds = [], units = [], volumes = [],
        inMM = [], notes = [];
    f.rows.forEach(function (row) {
      var name = row[cCol];
      var stock = C.cellNum(f, row, sCol);
      var desired = C.cellNum(f, row, dCol);
      var fixed = C.cellNum(f, row, vCol);
      var uL;
      if (fixed !== null) {
        uL = fixed;
      } else if (stock !== null && desired !== null) {
        if (!(stock > 0))
          throw new Error("'" + name + "' has a stock of " + stock +
                          '; it must be greater than 0');
        if (desired > stock)
          throw new Error("'" + name + "' wants " + desired +
                          ' from a stock of ' + stock +
                          '; it cannot be reached by dilution');
        uL = rxn * desired / stock;
      } else {
        throw new Error("'" + name + "' needs either a volume in '" + vCol +
                        "' or both '" + sCol + "' and '" + dCol + "'");
      }
      names.push(name);
      stocks.push(C.blank(C.cell(f, row, sCol)));
      desireds.push(C.blank(C.cell(f, row, dCol)));
      units.push(C.blank(C.cell(f, row, uCol)));
      volumes.push(C.pyRound(uL, 2));
      var keep = true;
      if (C.has(f, mCol)) {
        var k = row[mCol];
        keep = (k === null || k === undefined) ? true : Boolean(k);
      }
      inMM.push(keep);
      notes.push(C.blank(C.cell(f, row, nCol)));
    });

    var total = volumes.reduce(function (t, v) { return t + v; }, 0);
    var waterUL = C.pyRound(rxn - total, 2);
    if (waterUL < 0)
      throw new Error('The components come to ' + total.toFixed(2) +
                      ' uL, which is more than the ' + C.label(o, 'rxn_uL') +
                      ' uL reaction. Raise rxn_uL or lower a desired ' +
                      'concentration.');

    // Water first, the way it goes into the tube; total last
    names = [water].concat(names, ['Total']);
    stocks = [''].concat(stocks, ['']);
    desireds = [''].concat(desireds, ['']);
    units = [''].concat(units, ['']);
    volumes = [waterUL].concat(volumes, [C.pyRound(rxn, 2)]);
    inMM = [true].concat(inMM, [true]);
    notes = [''].concat(notes, ['']);

    var scale = samples * mmX;
    var mmUL = volumes.map(function (v, i) {
      return inMM[i] ? C.pyRound(v * scale, 2) : '-';
    });
    // The total is what goes in the master mix tube, scaled from the unrounded
    // sum so it is the volume to expect rather than the sum of rounded cells
    var kept = 0;
    for (var i = 0; i < volumes.length - 1; i++) if (inMM[i]) kept += volumes[i];
    mmUL[mmUL.length - 1] = C.pyRound(kept * scale, 2);

    var head = 'uL, MM (' + o.samples_label + ' x ' + o.mm_x_label + ')';
    var columns = ['Component', 'Stock', 'Desired', 'Unit', 'uL', head, 'Note'];
    var rows = names.map(function (n, i) {
      var r = {'Component': n, 'Stock': stocks[i], 'Desired': desireds[i],
               'Unit': units[i], 'uL': volumes[i], 'Note': notes[i]};
      r[head] = mmUL[i];
      return r;
    });
    if (!notes.some(function (n) { return typeof n === 'string' && n; })) {
      columns = columns.filter(function (c) { return c !== 'Note'; });
      rows.forEach(function (r) { delete r['Note']; });
    }
    return {columns: columns, rows: rows};
  };

  C['edms.bio.clone:rxn_mm'] = function (a) {
    return C._rxnMM(a.df, {
      rxn_uL: a.rxn_uL, samples: a.samples === undefined ? 1 : a.samples,
      mm_x: a.mm_x, water: a.water,
      samples_label: C.label(a, 'samples', 1),
      mm_x_label: C.label(a, 'mm_x', 1.1),
      component_col: a.component_col, stock_col: a.stock_col,
      desired_col: a.desired_col, unit_col: a.unit_col, uL_col: a.uL_col,
      mm_col: a.mm_col, note_col: a.note_col});
  };

  // ----------------------------------------------------------------- //
  // clone:anneal_mm
  // ----------------------------------------------------------------- //
  C['edms.bio.clone:anneal_mm'] = function (a) {
    var components = [
      {'Component': 'Top oligo (100 uM)', 'uL': C.g(a, 'top_uL', 1.0),
       'In MM': false},
      {'Component': 'Bottom oligo (100 uM)', 'uL': C.g(a, 'bottom_uL', 1.0),
       'In MM': false},
      {'Component': 'T4 DNA ligase buffer', 'Stock': C.g(a, 'buffer_stock', 10),
       'Desired': C.g(a, 'buffer_desired', 1), 'Unit': 'x'},
      {'Component': 'T4 PNK', 'uL': C.g(a, 'pnk_uL', 0.5), 'Note': 'add last'}
    ];
    var program = {columns: ['Step', 'Temperature', 'Time', 'Note'], rows: [
      {'Step': 1, 'Temperature': '37 C', 'Time': '30 min',
       'Note': 'phosphorylate'},
      {'Step': 2, 'Temperature': '95 C', 'Time': '5 min', 'Note': 'denature'},
      {'Step': 3, 'Temperature': '4 C', 'Time': 'hold',
       'Note': 'ramp at 0.2 C/sec so the duplexes anneal rather than snap back'}
    ]};
    return {__tables__: {
      'Anneal': C._rxnMM(components, {
        rxn_uL: a.rxn_uL === undefined ? 10 : a.rxn_uL,
        samples: a.oligos === undefined ? 1 : a.oligos,
        mm_x: a.mm_x === undefined ? 1.05 : a.mm_x,
        samples_label: C.label(a, 'oligos', 1),
        mm_x_label: C.label(a, 'mm_x', 1.05)}),
      'Anneal program': program}};
  };

  // ----------------------------------------------------------------- //
  // clone:golden_gate
  // ----------------------------------------------------------------- //
  C['edms.bio.clone:golden_gate'] = function (a) {
    var inserts = parseInt(a.inserts === undefined ? 3 : a.inserts, 10);
    var vectorStock = C.g(a, 'vector_stock', 100);
    var vectorNg = C.g(a, 'vector_ng', 100);
    var insertUL = C.g(a, 'insert_uL', null);
    var shared = a.shared_inserts === undefined ? false
      : Boolean(a.shared_inserts);
    var cycles = a.cycles === undefined ? 15 : a.cycles;
    if (!(inserts >= 1))
      throw new Error('inserts must be at least 1, got ' + inserts);
    if (!(vectorStock > 0))
      throw new Error('vector_stock must be greater than 0, got ' + vectorStock);

    var rows = [
      {'Component': 'FastDigest buffer', 'Stock': C.g(a, 'buffer_stock', 10),
       'Desired': C.g(a, 'buffer_desired', 1), 'Unit': 'x'},
      {'Component': 'ATP', 'Stock': C.g(a, 'atp_stock', 10),
       'Desired': C.g(a, 'atp_desired', 1), 'Unit': 'mM'}
    ];
    for (var i = 1; i <= inserts; i++) {
      var row = {'Component': 'Fragment ' + i, 'In MM': shared};
      if (insertUL) {  // falsy zero means the insert is a dilution instead
        row['uL'] = insertUL;
      } else {
        row['Stock'] = C.g(a, 'insert_stock', 10);
        row['Desired'] = C.g(a, 'insert_desired', 1);
        row['Unit'] = a.insert_unit === undefined ? 'x' : a.insert_unit;
      }
      rows.push(row);
    }
    rows.push({'Component': a.vector === undefined ? 'Vector' : a.vector,
               'uL': C.pyRound(vectorNg / vectorStock, 2),
               'Note': C.pyG(vectorNg) + ' ng at ' + C.pyG(vectorStock) +
                       ' ng/uL'});
    var enzyme = a.enzyme === undefined ? 'FastDigest Esp3I' : a.enzyme;
    rows.push({'Component': enzyme, 'Stock': C.g(a, 'enzyme_stock', 20),
               'Desired': C.g(a, 'enzyme_desired', 1), 'Unit': 'x'});
    rows.push({'Component': 'T4 DNA ligase', 'Stock': C.g(a, 'ligase_stock', 20),
               'Desired': C.g(a, 'ligase_desired', 1), 'Unit': 'x'});

    var program = {columns: ['Step', 'Temperature', 'Time', 'Note'], rows: [
      {'Step': 1, 'Temperature': '37 C', 'Time': '5 min', 'Note': 'digest'},
      {'Step': 2, 'Temperature': '16 C', 'Time': '5 min',
       'Note': 'ligate; repeat steps 1-2 for ' + C.label(a, 'cycles', 15) +
               ' cycles'},
      {'Step': 3, 'Temperature': '4 C', 'Time': 'hold', 'Note': ''}
    ]};
    return {__tables__: {
      'Golden Gate': C._rxnMM(rows, {
        rxn_uL: a.rxn_uL === undefined ? 10 : a.rxn_uL,
        samples: a.samples === undefined ? 1 : a.samples,
        mm_x: a.mm_x,
        samples_label: C.label(a, 'samples', 1),
        mm_x_label: C.label(a, 'mm_x', 1.1)}),
      'Golden Gate program': program}};
  };

  // ----------------------------------------------------------------- //
  // clone:gibson
  // ----------------------------------------------------------------- //
  C['edms.bio.clone:gibson'] = function (a) {
    var f = C.frame(a.df);
    var fragCol = a.fragment_col || 'Fragment';
    var kbCol = a.kb_col || 'kb';
    var concCol = a.conc_col || 'ng/uL';
    var pmolCol = a.pmol_col || 'pmol';
    var dilCol = a.dilution_col || 'Dilution';
    [fragCol, kbCol, concCol].forEach(function (col) {
      if (!C.has(f, col))
        throw new Error("gibson() needs a '" + col + "' column; df has " +
                        f.columns.join(', '));
    });
    if (!f.rows.length) throw new Error('gibson() needs at least one fragment');
    var rxn = C.g(a, 'rxn_uL', 15);
    var mmUL = C.g(a, 'mm_uL', null);
    if (mmUL === null || mmUL === undefined || isNaN(mmUL)) mmUL = rxn / 2;
    if (mmUL >= rxn)
      throw new Error('mm_uL (' + mmUL + ') must be less than rxn_uL (' +
                      rxn + ')');
    var pmol = C.g(a, 'pmol', 0.02);
    var pmolSmall = C.g(a, 'pmol_small', 0.05);
    var smallKb = C.g(a, 'small_kb', 1.0);
    var mmName = a.mm === undefined ? '2x NEBuilder HiFi MM' : a.mm;

    var names = [], kbs = [], pmols = [], ngs = [], concs = [], uLs = [];
    f.rows.forEach(function (row) {
      var name = row[fragCol];
      var kb = parseFloat(row[kbCol]);
      if (!(kb > 0))
        throw new Error("'" + name + "' has a length of " + kb +
                        ' kb; it must be greater than 0');
      var conc = C.cellNum(f, row, concCol);
      if (conc === null || !(conc > 0))
        throw new Error("'" + name + "' has no measured concentration. " +
                        'Nanodrop the fragments before planning the assembly; ' +
                        'the volumes cannot be guessed from the lengths alone.');
      var want = C.cellNum(f, row, pmolCol);
      if (want === null) want = kb < smallKb ? pmolSmall : pmol;
      var dilution = C.cellNum(f, row, dilCol);
      if (dilution === null) dilution = 1.0;
      if (!(dilution > 0))
        throw new Error("'" + name + "' has a dilution of " + dilution +
                        '; it must be greater than 0');
      var ng = want * kb * EDMS.consts.NG_PER_PMOL_KB;
      names.push(name);
      kbs.push(C.pyRound(kb, 3));
      pmols.push(C.pyRound(want, 4));
      ngs.push(C.pyRound(ng, 2));
      concs.push(C.pyRound(conc / dilution, 2));
      uLs.push(C.pyRound(ng / (conc / dilution), 2));
    });

    var used = uLs.reduce(function (t, v) { return t + v; }, 0);
    var waterUL = C.pyRound(rxn - mmUL - used, 2);
    if (waterUL < 0)
      throw new Error('The fragments come to ' + used.toFixed(2) +
                      ' uL, which leaves no room beside ' + C.label(a, 'mm_uL', mmUL) +
                      ' uL of master mix in a ' + C.label(a, 'rxn_uL', 15) +
                      ' uL reaction. Concentrate the fragments, lower the ' +
                      'pmol target, or scale the reaction up.');

    var rows = names.map(function (n, i) {
      return {'Fragment': n, 'kb': kbs[i], 'pmol': pmols[i], 'ng': ngs[i],
              'ng/uL': concs[i], 'uL': uLs[i]};
    });
    rows.push({'Fragment': mmName, 'kb': '', 'pmol': '', 'ng': '', 'ng/uL': '',
               'uL': C.pyRound(mmUL, 2)});
    rows.push({'Fragment': 'Nuclease-free water', 'kb': '', 'pmol': '', 'ng': '',
               'ng/uL': '', 'uL': waterUL});
    rows.push({'Fragment': 'Total', 'kb': '', 'pmol': '', 'ng': '', 'ng/uL': '',
               'uL': C.pyRound(rxn, 2)});
    return {columns: ['Fragment', 'kb', 'pmol', 'ng', 'ng/uL', 'uL'],
            rows: rows};
  };

  // ----------------------------------------------------------------- //
  // culture -- shared pieces
  // ----------------------------------------------------------------- //
  /* Python's format(x, ',.0f'): rounded to whole cells, grouped in threes.
     Only an error message needs it, but an error message about cell counts
     is unreadable without the commas. */
  C.pyComma = function (x, dp) {
    var s = Math.abs(C.pyRound(x, dp || 0)).toFixed(dp || 0);
    var dot = s.indexOf('.');
    var head = dot === -1 ? s : s.slice(0, dot);
    var tail = dot === -1 ? '' : s.slice(dot);
    head = head.replace(/\B(?=(\d{3})+(?!\d))/g, ',');
    return (x < 0 ? '-' : '') + head + tail;
  };

  /* A table of Component / Value / Unit, which is how every culture
     calculator reports itself. */
  C._cvu = function (rows) {
    return {columns: ['Component', 'Value', 'Unit'],
            rows: rows.map(function (r) {
              return {'Component': r[0], 'Value': r[1], 'Unit': r[2]};
            })};
  };

  C._confluentCells = function (vessel, confluency, recovery, override) {
    confluency = confluency === undefined ? 100 : confluency;
    recovery = (recovery === undefined || recovery === null) ? 1.0 : recovery;
    if (!(confluency > 0 && confluency <= 100))
      throw new Error('confluency must be a percent above 0 and at most 100, ' +
                      'got ' + confluency);
    if (!(recovery > 0 && recovery <= 1))
      throw new Error('recovery is a fraction above 0 and at most 1, got ' +
                      recovery);
    var name = C.canonical(vessel);
    var table = {};
    Object.keys(EDMS.consts.CONFLUENT_CELLS).forEach(function (k) {
      table[k] = EDMS.consts.CONFLUENT_CELLS[k];
    });
    Object.keys(override || {}).forEach(function (k) {
      table[C.canonical(k)] = override[k];
    });
    if (!Object.prototype.hasOwnProperty.call(table, name))
      throw new Error("No confluent cell count for '" + name +
                      "'. Known vessels: " + Object.keys(table).sort().join(', ') +
                      '. Pass confluent_cells_ to add it.');
    return table[name] * confluency / 100 * recovery;
  };

  C._perML = function (vessel, confluency, harvestML, override) {
    if (!(confluency > 0 && confluency <= 100))
      throw new Error('confluency must be a percent above 0 and at most 100, ' +
                      'got ' + confluency);
    if (harvestML !== undefined && harvestML !== null) {
      if (!(harvestML > 0))
        throw new Error('harvest_mL must be greater than 0, got ' + harvestML);
      return C._confluentCells(vessel, confluency) / harvestML;
    }
    var name = C.canonical(vessel);
    var table = {};
    Object.keys(EDMS.consts.CONFLUENT_PER_ML).forEach(function (k) {
      table[k] = EDMS.consts.CONFLUENT_PER_ML[k];
    });
    Object.keys(override || {}).forEach(function (k) {
      table[C.canonical(k)] = override[k];
    });
    if (!Object.prototype.hasOwnProperty.call(table, name))
      throw new Error("No confluent density for '" + name +
                      "'. Known vessels: " + Object.keys(table).sort().join(', ') +
                      '. Pass harvest_mL and it is worked out from the ' +
                      'confluent cell count instead, or pass confluent_per_mL, ' +
                      'or give the measured concentration per sample.');
    return table[name] * confluency / 100;
  };

  // ----------------------------------------------------------------- //
  // culture:drug_mm
  // ----------------------------------------------------------------- //
  C._drugMM = function (o) {
    var drugML = o.drug_media_mL, cellML = o.cell_media_mL;
    if (!(drugML > 0))
      throw new Error('drug_media_mL must be greater than 0, got ' + drugML);
    if (!(o.stock > 0))
      throw new Error('stock must be greater than 0, got ' + o.stock);
    if (o.desired > o.stock)
      throw new Error('desired (' + o.desired + ' ' + o.unit + ') is above ' +
                      'the stock concentration (' + o.stock + ' ' + o.unit +
                      '); it cannot be reached by dilution');
    var finalML = drugML + cellML;
    var concMM = o.desired * finalML / drugML;
    var totalML = drugML * o.wells * o.mm_x;
    var stockML = totalML * concMM / o.stock;
    return C._cvu([
      [o.drug + ', stock', o.stock_label, o.unit],
      [o.drug + ', desired', o.desired_label, o.unit],
      ['Drug media / well', o.drug_media_mL_label, 'mL'],
      ['Cell media / well', o.cell_media_mL_label, 'mL'],
      ['Wells', o.wells, 'wells'],
      ['MM', o.mm_x_label, 'x'],
      [o.drug + ', MM', C.pyRound(concMM, 4), o.unit],
      [o.drug + ' stock, MM', C.pyRound(stockML, 4), 'mL'],
      ['Media, MM', C.pyRound(totalML - stockML, 4), 'mL']]);
  };

  C['edms.bio.culture:drug_mm'] = function (a) {
    return C._drugMM({
      stock: C.g(a, 'stock', null), desired: C.g(a, 'desired', null),
      wells: parseInt(a.wells, 10),
      drug_media_mL: C.g(a, 'drug_media_mL', 1.0),
      cell_media_mL: C.g(a, 'cell_media_mL', 1.0),
      mm_x: C.g(a, 'mm_x', 1.1),
      drug: a.drug === undefined ? 'Drug' : a.drug,
      unit: a.unit === undefined ? 'ug/mL' : a.unit,
      stock_label: C.raw(a, 'stock'),
      desired_label: C.raw(a, 'desired'),
      drug_media_mL_label: C.raw(a, 'drug_media_mL', 1.0),
      cell_media_mL_label: C.raw(a, 'cell_media_mL', 1.0),
      mm_x_label: C.raw(a, 'mm_x', 1.1)});
  };

  // ----------------------------------------------------------------- //
  // culture:cells_mm
  // ----------------------------------------------------------------- //
  C._cellsMM = function (o) {
    if (!(o.flow_per_mL > 0))
      throw new Error('flow_per_mL must be greater than 0, got ' +
                      o.flow_per_mL);
    var cellsNeeded = o.cells_per_well * o.wells * o.mm_x;
    var spinML = cellsNeeded / o.flow_per_mL;
    var resuspendML = o.cell_media_mL * o.wells * o.mm_x;
    return C._cvu([
      ['Cells, flow', o.flow_per_mL_label, 'cells/mL'],
      ['Cells, desired', o.cells_per_well_label, 'cells/well'],
      ['Cell media / well', o.cell_media_mL_label, 'mL'],
      ['Wells', o.wells, 'wells'],
      ['MM', o.mm_x_label, 'x'],
      ['Cells needed, MM', C.pyRound(cellsNeeded, 1), 'cells'],
      ['Resuspension, MM', C.pyRound(cellsNeeded / resuspendML, 1), 'cells/mL'],
      ['Culture to spin, MM', C.pyRound(spinML, 4), 'mL'],
      ['Resuspend in, MM', C.pyRound(resuspendML, 4), 'mL']]);
  };

  // ----------------------------------------------------------------- //
  // culture:transduction
  // ----------------------------------------------------------------- //
  C['edms.bio.culture:transduction'] = function (a) {
    var wells = parseInt(a.wells, 10);
    var cellML = C.g(a, 'cell_media_mL', 1.0);
    var virusML = C.g(a, 'virus_mL_per_well', 1.0);
    // Polybrene is made up in the cell media but has to cover the virus
    // volume too, so the virus volume is what the cells are added on top of
    var polybrene = C._drugMM({
      stock: C.g(a, 'polybrene_stock', 10000),
      desired: C.g(a, 'polybrene_desired', 8),
      wells: wells, drug_media_mL: cellML, cell_media_mL: virusML,
      mm_x: C.g(a, 'mm_x_polybrene', 1.2), drug: 'Polybrene',
      unit: a.polybrene_unit === undefined ? 'ug/mL' : a.polybrene_unit,
      stock_label: C.raw(a, 'polybrene_stock', 10000),
      desired_label: C.raw(a, 'polybrene_desired', 8),
      drug_media_mL_label: C.raw(a, 'cell_media_mL', 1.0),
      cell_media_mL_label: C.raw(a, 'virus_mL_per_well', 1.0),
      mm_x_label: C.raw(a, 'mm_x_polybrene', 1.2)});
    var cells = C._cellsMM({
      flow_per_mL: C.g(a, 'flow_per_mL', null),
      cells_per_well: C.g(a, 'cells_per_well', null),
      wells: wells, cell_media_mL: cellML, mm_x: C.g(a, 'mm_x_cells', 1.1),
      flow_per_mL_label: C.raw(a, 'flow_per_mL'),
      cells_per_well_label: C.raw(a, 'cells_per_well'),
      cell_media_mL_label: C.raw(a, 'cell_media_mL', 1.0),
      mm_x_label: C.raw(a, 'mm_x_cells', 1.1)});
    return {__tables__: {'Polybrene': polybrene, 'Cells': cells}};
  };

  // ----------------------------------------------------------------- //
  // culture:sort_plan
  // ----------------------------------------------------------------- //
  C['edms.bio.culture:sort_plan'] = function (a) {
    var target = parseInt(a.target_clones, 10);
    var eff = C.g(a, 'plating_efficiency', 0.1);
    var fraction = C.g(a, 'conditioned_fraction', 0.5);
    if (!(target >= 1))
      throw new Error('target_clones must be at least 1, got ' + target);
    if (!(eff > 0 && eff <= 1))
      throw new Error('plating_efficiency is a fraction above 0 and at most ' +
                      '1, got ' + eff);
    if (!(fraction >= 0 && fraction <= 1))
      throw new Error('conditioned_fraction is a fraction between 0 and 1, ' +
                      'got ' + fraction);
    var s = C.spec(a.plate_format === undefined ? '96-well' : a.plate_format);
    var mediaML = (a.media_mL === undefined || a.media_mL === null ||
                   a.media_mL === '') ? s.media_mL : parseFloat(a.media_mL);
    var mediaLabel = (a.media_mL === undefined || a.media_mL === null ||
                      a.media_mL === '') ? s.media_mL : C.raw(a, 'media_mL');
    var wellsNeeded = Math.ceil(target / eff);
    var plateCount = Math.ceil(wellsNeeded / s.wells);
    var wellsSorted = plateCount * s.wells;
    var totalML = wellsSorted * mediaML;
    return C._cvu([
      ['Clones wanted', target, 'clones'],
      ['Plating efficiency', C.raw(a, 'plating_efficiency', 0.1),
       'fraction'],
      ['Format', s.plate_format, '-'],
      ['Wells needed', wellsNeeded, 'wells'],
      ['Plates to sort', plateCount, 'plates'],
      ['Wells sorted', wellsSorted, 'wells'],
      ['Clones expected', C.pyRound(wellsSorted * eff, 1), 'clones'],
      ['Media per well', mediaLabel, 'mL'],
      ['Media total', C.pyRound(totalML, 1), 'mL'],
      ['Conditioned media', C.pyRound(totalML * fraction, 1), 'mL'],
      ['Fresh media', C.pyRound(totalML * (1 - fraction), 1), 'mL']]);
  };

  // ----------------------------------------------------------------- //
  // culture:freeze_down
  // ----------------------------------------------------------------- //
  C['edms.bio.culture:freeze_down'] = function (a) {
    var harvestML = C.g(a, 'harvest_mL', null);
    var perVial = C.g(a, 'cells_per_vial', 1e6);
    var dmso = C.g(a, 'dmso_pct', 10);
    var recovery = C.g(a, 'recovery', 1.0);
    if (!(harvestML > 0))
      throw new Error('harvest_mL must be greater than 0, got ' + harvestML);
    if (!(perVial > 0))
      throw new Error('cells_per_vial must be greater than 0, got ' + perVial);
    if (!(dmso >= 0 && dmso < 100))
      throw new Error('dmso_pct must be at least 0 and below 100, got ' + dmso);

    var measured = C.g(a, 'cells_per_mL', null);
    var conc, total, source;
    if (measured !== null && !isNaN(measured)) {
      if (!(measured > 0))
        throw new Error('cells_per_mL must be greater than 0, got ' + measured);
      conc = measured; source = 'measured'; total = conc * harvestML;
    } else if (a.vessel !== undefined && a.vessel !== null && a.vessel !== '' &&
               a.confluency !== undefined && a.confluency !== null &&
               a.confluency !== '') {
      var confluency = parseFloat(a.confluency);
      total = C._confluentCells(a.vessel, confluency, recovery,
                                a.confluent_cells_);
      conc = total / harvestML;
      source = C.pyG(confluency) + '% of a ' + C.canonical(a.vessel);
      if (recovery !== 1.0) source += ' at ' + C.pyG(recovery) + ' recovery';
    } else {
      throw new Error('freeze_down() needs either cells_per_mL, or both ' +
                      'vessel and confluency to estimate one');
    }

    var vials = Math.floor(total / perVial);
    if (vials < 1)
      throw new Error(C.pyComma(total, 0) + ' cells is less than the ' +
                      C.pyComma(perVial, 0) + ' wanted in one vial. Grow the ' +
                      'clone further, or lower cells_per_vial.');
    var mediumML = vials * C.g(a, 'vial_mL', 1.0);
    return C._cvu([
      ['Cells, source', source, '-'],
      ['Cells', C.pyRound(conc, 1), 'cells/mL'],
      ['Harvest', C.raw(a, 'harvest_mL'), 'mL'],
      ['Cells, total', C.pyRound(total, 1), 'cells'],
      ['Cells / vial', C.raw(a, 'cells_per_vial', 1e6), 'cells'],
      ['Vials', vials, 'vials'],
      ['Cells left over', C.pyRound(total - vials * perVial, 1), 'cells'],
      ['Freezing medium', C.pyRound(mediumML, 2), 'mL'],
      ['DMSO', C.pyRound(mediumML * dmso / 100, 2), 'mL'],
      ['Media / FBS', C.pyRound(mediumML * (1 - dmso / 100), 2), 'mL']]);
  };

  // ----------------------------------------------------------------- //
  // culture:seed_plan -- one arm of an experiment, per sample
  // ----------------------------------------------------------------- //
  C.merge = function (a, b) {
    var out = {};
    Object.keys(a || {}).forEach(function (k) { out[k] = a[k]; });
    Object.keys(b || {}).forEach(function (k) { out[k] = b[k]; });
    return out;
  };

  C._seedPlan = function (f, o) {
    var confCol = o.confluency_col || 'Confluency %';
    var vesCol = o.vessel_col || 'Vessel';
    var perCol = o.per_mL_col || 'Cells, flow';
    [confCol, vesCol].forEach(function (col) {
      if (!C.has(f, col))
        throw new Error("seed_plan() needs a '" + col + "' column; df has " +
                        f.columns.join(', '));
    });
    if (!(o.wells > 0))
      throw new Error('wells must be greater than 0, got ' + o.wells);

    var cellsNeeded = o.cells_per_well * o.wells * o.mm_x;
    var resuspendML = o.cell_media_mL * o.wells * o.mm_x;
    var spinCol = o.arm + ' spin, mL';
    var resCol = o.arm + ' resuspend, mL';

    var rows = f.rows.map(function (row) {
      var out = {};
      f.columns.forEach(function (c) { out[c] = row[c]; });
      var measured = C.cellNum(f, row, perCol);
      var conc;
      if (measured !== null) {
        if (!(measured > 0))
          throw new Error("'" + perCol + "' must be greater than 0, got " +
                          measured);
        conc = measured;
      } else {
        conc = C._perML(row[vesCol], parseFloat(row[confCol]), null,
                        o.confluent_per_mL);
      }
      var scale = 1.0;
      if (o.scale_col) {
        var sv = C.cellNum(f, row, o.scale_col);
        if (sv !== null) scale = sv;
      }
      out[spinCol] = C.pyRound(cellsNeeded * scale / conc, 4);
      out[resCol] = C.pyRound(resuspendML * scale, 4);
      return out;
    });

    var columns = f.columns.concat([spinCol, resCol]);
    /* Python selects the spin columns by suffix, which matches 'Sum spin, mL'
       as well. It cannot bite at two arms, which is all competitive_growth
       asks for, and a port that quietly disagreed with the package would be
       worse than one that follows it here. */
    var spinCols = columns.filter(function (c) {
      return c.length > 9 && c.slice(-9) === ' spin, mL';
    });
    if (spinCols.length > 1) {
      rows.forEach(function (r) {
        var t = 0;
        spinCols.forEach(function (c) {
          var v = C.num(r[c]);
          if (v !== null) t += v;
        });
        r['Sum spin, mL'] = C.pyRound(t, 4);
      });
      // Keep the total last however many arms have been added
      columns = columns.filter(function (c) { return c !== 'Sum spin, mL'; });
      columns.push('Sum spin, mL');
    }
    return {columns: columns, rows: rows};
  };

  // ----------------------------------------------------------------- //
  // culture:competitive_growth
  // ----------------------------------------------------------------- //
  C['edms.bio.culture:competitive_growth'] = function (a) {
    var f = C.frame(a.df);
    if (!f.rows.length)
      throw new Error('competitive_growth() needs at least one sample');
    var cgCells = C.g(a, 'cg_cells_per_well', 1.2e5);
    var cgWells = parseInt(a.cg_wells === undefined ? 3 : a.cg_wells, 10);
    var tCells = C.g(a, 'transduce_cells_per_well', 1.8e5);
    var tWells = parseInt(a.transduce_wells === undefined ? 3
                          : a.transduce_wells, 10);
    var cellML = C.g(a, 'cell_media_mL', 0.25);
    var mmX = C.g(a, 'mm_x', 1.2);
    var compMMX = C.g(a, 'competitor_mm_x', 1.2);
    var extra = parseInt(a.extra_samples === undefined ? 1 : a.extra_samples, 10);
    var competitor = a.competitor === undefined ? 'LCG' : a.competitor;
    var vessel = a.competitor_vessel === undefined ? '15-cm'
      : a.competitor_vessel;
    var confluency = C.g(a, 'competitor_confluency', 100);

    var shared = {cell_media_mL: cellML, mm_x: mmX,
                  confluency_col: a.confluency_col, vessel_col: a.vessel_col,
                  per_mL_col: a.per_mL_col, scale_col: a.scale_col,
                  confluent_per_mL: a.confluent_per_mL};
    var plan = f;
    if (tWells)
      plan = C._seedPlan(plan, C.merge(shared, {cells_per_well: tCells,
                                                wells: tWells,
                                                arm: 'Transduce'}));
    plan = C._seedPlan(plan, C.merge(shared, {cells_per_well: cgCells,
                                              wells: cgWells, arm: 'CG'}));

    // The competitor is made up once for every sample, plus the day 0 mix
    var samples = f.rows.length + extra;
    var cellsNeeded = cgCells * cgWells * mmX * samples * compMMX;
    var resuspendML = cellML * cgWells * mmX * samples * compMMX;
    var measured = C.g(a, 'competitor_per_mL', null);
    var conc, source;
    if (measured !== null && !isNaN(measured)) {
      if (!(measured > 0))
        throw new Error('competitor_per_mL must be greater than 0, got ' +
                        measured);
      conc = measured;
      source = 'measured';
    } else {
      conc = C._perML(vessel, confluency, null, a.confluent_per_mL);
      source = C.pyG(confluency) + '% of a ' + vessel;
    }
    var competitorTable = C._cvu([
      [competitor + ', source', source, '-'],
      [competitor + ', cells', C.pyRound(conc, 1), 'cells/mL'],
      ['Cells, desired', C.raw(a, 'cg_cells_per_well', 1.2e5), 'cells/well'],
      ['Wells / sample', C.raw(a, 'cg_wells', 3), 'wells'],
      ['Samples', samples, 'samples'],
      ['MM', C.pyRound(mmX * compMMX, 4), 'x'],
      [competitor + ' needed, MM', C.pyRound(cellsNeeded, 1), 'cells'],
      [competitor + ' to spin, MM', C.pyRound(cellsNeeded / conc, 4), 'mL'],
      ['Resuspend in, MM', C.pyRound(resuspendML, 4), 'mL']]);

    var tables = {'Seed plan': plan};
    tables[competitor] = competitorTable;
    return {__tables__: tables};
  };

  // ----------------------------------------------------------------- //
  // ngs:coverage
  // ----------------------------------------------------------------- //
  C['edms.bio.ngs:coverage'] = function (a) {
    var genotypes = C.g(a, 'genotypes', null);
    var samples = parseInt(a.samples === undefined ? 1 : a.samples, 10);
    var cellCoverage = C.g(a, 'cell_coverage', 1000);
    var ugPerCell = C.g(a, 'ug_gDNA_per_cell', 6 * Math.pow(10, -6));
    var ploidy = C.g(a, 'ploidy_per_cell', 2);
    var umiCoverage = C.g(a, 'umi_coverage', 5);
    var ngPerPcr = C.g(a, 'ng_per_pcr', 200);
    if (!(ngPerPcr > 0))
      throw new Error('ng_per_pcr must be greater than 0, got ' + ngPerPcr);

    var ug = genotypes * cellCoverage * ugPerCell;
    var molecules = genotypes * cellCoverage * ploidy;
    var reads = molecules * umiCoverage;
    var readsAll = reads * samples;
    var reactions = Math.ceil(ug * 1000 / ngPerPcr);

    var rows = [
      ['Genotypes', C.raw(a, 'genotypes'), 'genotypes'],
      ['Cells per genotype', C.raw(a, 'cell_coverage', 1000), 'cells'],
      ['Samples', C.raw(a, 'samples', 1), 'samples'],
      ['gDNA per sample', C.pyRound(ug, 3), 'ug'],
      ['gDNA per sample', C.pyRound(ug * 1000, 1), 'ng'],
      ['PCR1 reactions per sample', reactions,
       'reactions at ' + C.label(a, 'ng_per_pcr', 200) + ' ng'],
      ['Molecules per sample', C.pyRound(molecules / 1e6, 2), 'M'],
      ['Reads per sample', C.pyRound(reads / 1e6, 2), 'M'],
      ['Reads for all samples', C.pyRound(readsAll / 1e6, 2), 'M']];
    return {columns: ['Quantity', 'Value', 'Unit'],
            rows: rows.map(function (r) {
              return {'Quantity': r[0], 'Value': r[1], 'Unit': r[2]};
            })};
  };

  // ----------------------------------------------------------------- //
  // plates:PE3  (wraps edms.bio.transfect:PE3)
  // ----------------------------------------------------------------- //
  /* pandas merges two frames on one key by keeping the left frame's columns
     in order, then the right frame's minus the key, suffixing any name that
     appears in both. A nanodrop sheet with extra columns hits that path, so
     it is mirrored rather than assumed away. */
  C._mergeOn = function (left, right, key) {
    var shared = left.columns.filter(function (c) {
      return c !== key && right.columns.indexOf(c) !== -1;
    });
    var cols = [], map = [];
    left.columns.forEach(function (c) {
      var name = shared.indexOf(c) === -1 ? c : c + '_x';
      cols.push(name); map.push(['l', c, name]);
    });
    right.columns.forEach(function (c) {
      if (c === key) return;
      var name = shared.indexOf(c) === -1 ? c : c + '_y';
      cols.push(name); map.push(['r', c, name]);
    });
    var rows = [];
    left.rows.forEach(function (lr) {
      right.rows.forEach(function (rr) {
        if (lr[key] !== rr[key]) return;
        var out = {};
        map.forEach(function (m) {
          out[m[2]] = (m[0] === 'l' ? lr : rr)[m[1]];
        });
        rows.push(out);
      });
    });
    return {columns: cols, rows: rows};
  };

  C._rename = function (f, mapping) {
    var cols = f.columns.map(function (c) {
      return Object.prototype.hasOwnProperty.call(mapping, c) ? mapping[c] : c;
    });
    var rows = f.rows.map(function (r) {
      var out = {};
      f.columns.forEach(function (c, i) { out[cols[i]] = r[c]; });
      return out;
    });
    return {columns: cols, rows: rows};
  };

  C._transfectPE3 = function (a) {
    var plasmids = C.frame(a.plasmids);
    var epegRNAs = C.frame(a.epegRNAs);
    var ngRNAs = C.frame(a.ngRNAs);
    var numCol = a.pegRNA_number_col || 'pegRNA_number';
    var eNameCol = a.epegRNAs_name_col || 'Name';
    var nNameCol = a.ngRNAs_name_col || 'Name';
    var pCol = a.plasmid_col || 'Plasmid';
    var dCol = a.description_col || 'Description';
    var cCol = a.colony_col || 'Colony';
    var concCol = a.ng_uL_col || 'ng/uL';
    var PE = a.PE_plasmid === undefined ? 'pMUZ86.7' : a.PE_plasmid;
    var reps = parseInt(a.reps === undefined ? 3 : a.reps, 10);
    var mmX = C.g(a, 'mm_x', 1.1);

    var pick = function (name) {
      return {columns: plasmids.columns,
              rows: plasmids.rows.filter(function (r) {
                return r[pCol] === name;
              })};
    };

    var tubeA = null;
    epegRNAs.rows.forEach(function (er) {
      var number = er[numCol], epegRNA = er[eNameCol];
      ngRNAs.rows.forEach(function (nr) {
        if (nr[numCol] !== number) return;
        var ngRNA = nr[nNameCol];
        var eMap = {}, nMap = {};
        eMap[pCol] = 'epegRNA';
        eMap[cCol] = 'epegRNA ' + cCol;
        eMap[dCol] = 'epegRNA ' + dCol;
        eMap[concCol] = 'epegRNA ' + concCol;
        nMap[pCol] = 'ngRNA';
        nMap[cCol] = 'ngRNA ' + cCol;
        nMap[dCol] = 'ngRNA ' + dCol;
        nMap['ng/uL'] = 'ngRNA ' + concCol;   // the Python hard-codes this one
        var ef = C._rename(pick(epegRNA), eMap);
        var nf = C._rename(pick(ngRNA), nMap);
        ef = {columns: ef.columns.concat(['ngRNA']),
              rows: ef.rows.map(function (r) {
                var o = C.merge(r, {}); o['ngRNA'] = ngRNA; return o;
              })};
        var merged = C._mergeOn(ef, nf, 'ngRNA');
        if (tubeA === null) tubeA = merged;
        else tubeA = {columns: tubeA.columns,
                      rows: tubeA.rows.concat(merged.rows)};
      });
    });
    if (tubeA === null || !tubeA.rows.length)
      throw new Error('No epegRNA and ngRNA pair shares a ' + numCol +
                      ', so there is nothing to transfect.');

    var peRows = plasmids.rows.filter(function (r) { return r[pCol] === PE; });
    if (!peRows.length)
      throw new Error("No plasmid named '" + PE + "' in the plasmids table.");
    var peColony = peRows[0][cCol], peConc = peRows[0][concCol];

    var epegNg = C.g(a, 'epegRNA_ng', 66);
    var ngNg = C.g(a, 'ngRNA_ng', 22);
    var peNg = C.g(a, 'PE_ng', 200);
    var wellUL = C.g(a, 'well_uL', 10);
    var half = C.pyRound(mmX * reps * wellUL / 2, 2);

    var cols = tubeA.columns.concat(
      ['PE', 'PE ' + cCol, 'PE ' + concCol, 'epegRNA uL', 'ngRNA uL', 'PE uL',
       'Optimem uL', 'L2000 uL', 'Total uL', 'Tube']);
    var body = tubeA.rows.map(function (r, i) {
      var o = C.merge(r, {});
      o['PE'] = PE;
      o['PE ' + cCol] = peColony;
      o['PE ' + concCol] = peConc;
      o['epegRNA uL'] = C.pyRound(mmX * reps * epegNg / C.num(r['epegRNA ' + concCol]), 2);
      o['ngRNA uL'] = C.pyRound(mmX * reps * ngNg / C.num(r['ngRNA ' + concCol]), 2);
      o['PE uL'] = C.pyRound(mmX * reps * peNg / C.num(peConc), 2);
      // The Optimem is what is left of the half reaction once the three
      // plasmids are in, worked out from the rounded volumes the way pandas
      // does it column by column
      o['Optimem uL'] = C.pyRound(mmX * reps * wellUL / 2 - o['epegRNA uL'] -
                                  o['ngRNA uL'] - o['PE uL'], 2);
      o['L2000 uL'] = 0;
      o['Total uL'] = half;
      o['Tube'] = 'A' + (i + 1);
      return o;
    });

    var n = body.length;
    var tubeB = {};
    cols.forEach(function (c) { tubeB[c] = null; });
    /* Straight from edms.bio.transfect:PE3. Its Optimem is 9/10 of
       mm_x*reps*mm_x*n rather than of the half reaction, so it does not add
       up to Total with the L2000 beside it. Mirrored rather than corrected:
       the page has to agree with the package, and a number that differs
       between the two is the one thing this file must never produce. */
    tubeB['Optimem uL'] = C.pyRound(9 / 10 * mmX * reps * mmX * n, 2);
    tubeB['L2000 uL'] = C.pyRound(1 / 10 * wellUL / 2 * mmX * reps * mmX * n, 2);
    tubeB['Total uL'] = C.pyRound(wellUL / 2 * mmX * reps * mmX * n, 2);
    tubeB['Tube'] = 'B';
    // concat makes the all-integer L2000 column a float once tube B joins it
    body.forEach(function (r) { r['L2000 uL'] = 0.0; });
    var transfection = {columns: cols, rows: body.concat([tubeB])};

    // Where each Tube A sits on the plate, skipping the outer wells
    var rowNames = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'];
    var placed = [];
    var plateI = 1, rowI = 1, colI = 1;
    for (var i = 0; i < n; i++) {
      if (colI >= 12 - 1) {
        if (rowI >= rowNames.length - 1 - reps) {
          rowI = 1; colI = 1; plateI += 1;
        } else {
          rowI += reps; colI = 1;
        }
      }
      placed.push({plate: plateI, row: rowNames[rowI], column: colI + 1,
                   Tube: body[i]['Tube']});
      colI += 1;
    }
    return {transfection: transfection, placed: placed};
  };

  /* pivot_table(values='Tube', index=['plate','row'], columns='column'):
     one row per (plate, row) in sorted order, one column per column number.
     interactive.table_data() carries a frame's columns and not its index, so
     the plate and row labels are not in the payload and are not here either. */
  C._pivotPlate = function (placed) {
    var colNums = [];
    placed.forEach(function (p) {
      if (colNums.indexOf(p.column) === -1) colNums.push(p.column);
    });
    colNums.sort(function (x, y) { return x - y; });
    var keys = [], seen = {};
    placed.forEach(function (p) {
      var k = p.plate + '\u0000' + p.row;
      if (!seen[k]) { seen[k] = true; keys.push([p.plate, p.row, k]); }
    });
    keys.sort(function (x, y) {
      if (x[0] !== y[0]) return x[0] - y[0];
      return x[1] < y[1] ? -1 : (x[1] > y[1] ? 1 : 0);
    });
    var rows = keys.map(function (k) {
      var out = {};
      colNums.forEach(function (c) {
        var hit = null;
        placed.forEach(function (p) {
          if (hit === null && p.plate === k[0] && p.row === k[1] &&
              p.column === c) hit = p.Tube;
        });
        out[String(c)] = hit;
      });
      return out;
    });
    return {columns: colNums.map(String), rows: rows};
  };

  C['edms.work.plates:PE3'] = function (a) {
    var name = C.canonical(a.plate_format === undefined
                           ? EDMS.consts.PE3_ANCHOR : a.plate_format);
    var f = C.scale(name, EDMS.consts.PE3_ANCHOR);
    var scaled = C.merge(a, {
      epegRNA_ng: C.pyRound(C.g(a, 'epegRNA_ng', 66) * f, 2),
      ngRNA_ng: C.pyRound(C.g(a, 'ngRNA_ng', 22) * f, 2),
      PE_ng: C.pyRound(C.g(a, 'PE_ng', 200) * f, 2),
      well_uL: C.pyRound(C.g(a, 'well_uL', 10) * f, 2)});
    var out = C._transfectPE3(scaled);
    var tables = {'Transfection': out.transfection};
    if (name === EDMS.consts.PE3_ANCHOR) {
      tables['96-well plates'] = C._pivotPlate(out.placed);
    } else {
      /* transfect.PE3 lays conditions out on a fixed 8x12 grid with the outer
         wells empty. That map describes a 96-well plate and nothing else. */
      tables[EDMS.consts.LAYOUT_NOTE + ' (' + name + ')'] = {
        columns: ['Note'],
        rows: [{'Note': 'No plate layout: edms.bio.transfect:PE3 lays out an ' +
                        '8x12 grid, which does not describe a ' + name +
                        ' plate. The volumes above are scaled correctly; ' +
                        'assign wells yourself.'}]};
    }
    return {__tables__: tables};
  };

  // ----------------------------------------------------------------- //
  // ngs:pcrs -- plate layouts, master mixes and thermocycler programs
  // ----------------------------------------------------------------- //
  /* Python's repr() of a tuple. edms.work.plan names a step's tables after
     the keys of whatever the calculator returned, and pcr_mm() keys its
     master mixes by a (FWD, REV) tuple, so the table is called
     "('P5-FWD', 'P7-REV')" -- quotes, comma, space and all. */
  C.pyRepr = function (v) {
    if (typeof v === 'string') {
      var q = (v.indexOf("'") !== -1 && v.indexOf('"') === -1) ? '"' : "'";
      var out = v.replace(/\\/g, '\\\\');
      if (q === "'") out = out.replace(/'/g, "\\'");
      return q + out + q;
    }
    if (v === null || v === undefined) return 'None';
    if (v === true) return 'True';
    if (v === false) return 'False';
    if (typeof v === 'number') return Number.isInteger(v) ? String(v)
                                                          : C.pyFloat(v);
    return String(v);
  };

  C.pyTuple = function (vals) {
    return '(' + vals.map(C.pyRepr).join(', ') + (vals.length === 1 ? ',' : '') +
      ')';
  };

  /* sum([a, b, c] * n) in Python is the list repeated n times and added left
     to right, not (a + b + c) * n. The two are different doubles, and the
     difference survives rounding to two decimals often enough to matter. */
  C._sumRep = function (parts, times) {
    var t = 0;
    for (var i = 0; i < times; i++)
      for (var j = 0; j < parts.length; j++) t += parts[j];
    return t;
  };

  C._uniqueTuples = function (f, cols) {
    var seen = {}, out = [];
    f.rows.forEach(function (r) {
      var vals = cols.map(function (c) { return r[c]; });
      var k = JSON.stringify(vals);
      if (!seen[k]) { seen[k] = true; out.push(vals); }
    });
    return out;
  };

  C._vcsOrdered = function (f, cols) {
    var counts = {};
    f.rows.forEach(function (r) {
      var k = JSON.stringify(cols.map(function (c) { return r[c]; }));
      counts[k] = (counts[k] || 0) + 1;
    });
    return C._uniqueTuples(f, cols).map(function (vals) {
      return {key: vals, count: counts[JSON.stringify(vals)]};
    });
  };

  C._groupBoundaries = function (nums) {
    if (!nums.length) return [];
    var sorted = nums.slice().sort(function (a, b) { return a - b; })
      .filter(function (v, i, arr) { return i === 0 || v !== arr[i - 1]; });
    var groups = [], start = sorted[0];
    for (var i = 1; i < sorted.length; i++) {
      if (sorted[i] - sorted[i - 1] > 1) {
        groups.push([start, sorted[i - 1]]);
        start = sorted[i];
      }
    }
    groups.push([start, sorted[sorted.length - 1]]);
    return groups;
  };

  C._minSec = function (decimalMinutes) {
    var minutes = Math.trunc(decimalMinutes);
    var seconds = Math.round(C.pyRound((decimalMinutes - minutes) * 60, 0));
    return [minutes, seconds];
  };

  var MM_COLS = ['Component', 'Stock', 'Desired', 'Unit', 'uL', 'uL MM'];

  C._addSyber = function (mm, template, primersN, totalUL, mmX, stock,
                          desired) {
    var syberUL = desired / stock * totalUL;
    var rows = mm.rows.map(function (r) { return C.merge(r, {}); });
    var h2o = null;
    rows.forEach(function (r) {
      if (h2o === null && r['Component'] === 'Nuclease-free H2O') h2o = r;
    });
    if (h2o) {
      h2o['uL'] = C.pyRound(h2o['uL'] - syberUL, 2);
      h2o['uL MM'] = C.pyRound(h2o['uL MM'] - syberUL * primersN * mmX, 2);
      if (h2o['uL'] < 0)
        throw new Error('Not enough Nuclease-free H2O for ' + C.label({}, 'x', stock) +
                        'x SYBR Green DNA Stain (' + C.pyRound(syberUL, 2) +
                        ' uL); decrease template or primer volumes.');
    }
    var at = -1;
    rows.forEach(function (r, i) {
      if (at === -1 && r['Component'] === template) at = i;
    });
    var row = {'Component': stock + 'x SYBR Green DNA Stain', 'Stock': stock,
               'Desired': desired, 'Unit': 'x', 'uL': C.pyRound(syberUL, 2),
               'uL MM': C.pyRound(syberUL * primersN * mmX, 2)};
    rows.splice(at === -1 ? rows.length : at, 0, row);
    return {columns: MM_COLS, rows: rows};
  };

  C._pcrMM = function (primers, o) {
    var out = [];
    primers.forEach(function (pair) {
      var fwd = pair.key[0], rev = pair.key[1], n = pair.count;
      var parts = [o.Q5_mm_x_desired / o.Q5_mm_x_stock,
                   o.dNTP_mM_desired / o.dNTP_mM_stock,
                   o.fwd_uM_desired / o.fwd_uM_stock,
                   o.rev_uM_desired / o.rev_uM_stock,
                   o.template_uL / o.total_uL,
                   o.Q5_U_uL_desired / o.Q5_U_uL_stock];
      var h2o = o.total_uL - C._sumRep(parts, o.total_uL);
      var uL = [h2o,
                o.Q5_mm_x_desired / o.Q5_mm_x_stock * o.total_uL,
                o.dNTP_mM_desired / o.dNTP_mM_stock * o.total_uL,
                o.fwd_uM_desired / o.fwd_uM_stock * o.total_uL,
                o.rev_uM_desired / o.rev_uM_stock * o.total_uL,
                o.template_uL,
                o.Q5_U_uL_desired / o.Q5_U_uL_stock * o.total_uL,
                o.total_uL];
      var names = ['Nuclease-free H2O', o.Q5_mm_x_stock + 'x Q5 Reaction Buffer',
                   'dNTPs', fwd, rev, o.template, 'Q5 Polymerase', 'Total'];
      var stocks = ['', o.Q5_mm_x_stock, o.dNTP_mM_stock, o.fwd_uM_stock,
                    o.rev_uM_stock, '', o.Q5_U_uL_stock, ''];
      var desireds = ['', o.Q5_mm_x_desired, o.dNTP_mM_desired,
                      o.fwd_uM_desired, o.rev_uM_desired, '',
                      o.Q5_U_uL_desired, ''];
      var units = ['', 'x', 'mM', 'uM', 'uM', '', 'U/uL', ''];
      var rows = names.map(function (name, i) {
        return {'Component': name, 'Stock': stocks[i], 'Desired': desireds[i],
                'Unit': units[i], 'uL': C.pyRound(uL[i], 2),
                'uL MM': C.pyRound(uL[i] * n * o.mm_x, 2)};
      });
      var table = {columns: MM_COLS, rows: rows};
      if (o.syber)
        table = C._addSyber(table, o.template, n, o.total_uL, o.mm_x,
                            o.syber_x_stock, o.syber_x_desired);
      out.push({key: [fwd, rev], table: table});
    });
    return out;
  };

  C._pcrMMUltra = function (primers, o) {
    var out = [];
    primers.forEach(function (pair) {
      var fwd = pair.key[0], rev = pair.key[1], n = pair.count;
      var parts = [o.Q5_mm_x_desired / o.Q5_mm_x_stock,
                   o.fwd_uM_desired / o.fwd_uM_stock,
                   o.rev_uM_desired / o.rev_uM_stock,
                   o.template_uL / o.total_uL];
      var h2o = o.total_uL - C._sumRep(parts, o.total_uL);
      var uL = [h2o,
                o.Q5_mm_x_desired / o.Q5_mm_x_stock * o.total_uL,
                o.fwd_uM_desired / o.fwd_uM_stock * o.total_uL,
                o.rev_uM_desired / o.rev_uM_stock * o.total_uL,
                o.template_uL,
                o.total_uL];
      var names = ['Nuclease-free H2O',
                   'NEBNext Ultra II Q5 ' + o.Q5_mm_x_stock + 'x MM',
                   fwd, rev, o.template, 'Total'];
      var stocks = ['', o.Q5_mm_x_stock, o.fwd_uM_stock, o.rev_uM_stock, '', ''];
      var desireds = ['', o.Q5_mm_x_desired, o.fwd_uM_desired,
                      o.rev_uM_desired, '', ''];
      var units = ['', 'x', 'uM', 'uM', '', ''];
      var rows = names.map(function (name, i) {
        return {'Component': name, 'Stock': stocks[i], 'Desired': desireds[i],
                'Unit': units[i], 'uL': C.pyRound(uL[i], 2),
                'uL MM': C.pyRound(uL[i] * n * o.mm_x, 2)};
      });
      var table = {columns: MM_COLS, rows: rows};
      if (o.syber)
        table = C._addSyber(table, o.template, n, o.total_uL, o.mm_x,
                            o.syber_x_stock, o.syber_x_desired);
      out.push({key: [fwd, rev], table: table});
    });
    return out;
  };

  C._thermocycler = function (f, o) {
    var n = String(o.n);
    var cycles = o.cycles;
    var fwdCol = o.pcr_fwd_col, revCol = o.pcr_rev_col;
    if (!o.umi) {
      if (n === '1') {
        if (cycles === null || cycles === undefined || cycles === '') cycles = 30;
        fwdCol = fwdCol || 'PCR1 FWD'; revCol = revCol || 'PCR1 REV';
      } else if (n === '2') {
        if (cycles === null || cycles === undefined || cycles === '') cycles = 8;
        fwdCol = fwdCol || 'PCR2 FWD'; revCol = revCol || 'PCR2 REV';
        C._setCol(f, fwdCol, 'PCR2-FWD');
      } else throw new Error('n must be 1 or 2');
    } else {
      if (n === '1') {
        if (cycles === null || cycles === undefined || cycles === '') cycles = 3;
        fwdCol = fwdCol || 'PCR1 FWD'; revCol = revCol || 'PCR1 REV';
      } else if (n === '1.5') {
        if (cycles === null || cycles === undefined || cycles === '') cycles = 30;
        fwdCol = fwdCol || 'PCR1 FWD'; revCol = revCol || 'PCR1 REV';
        C._setCol(f, fwdCol, 'P5-FWD');
        C._setCol(f, revCol, 'P7-REV');
      } else if (n === '2') {
        if (cycles === null || cycles === undefined || cycles === '') cycles = 8;
        fwdCol = fwdCol || 'PCR2 FWD'; revCol = revCol || 'PCR2 REV';
        C._setCol(f, fwdCol, 'PCR2-FWD');
      } else throw new Error('n must be 1, 1.5 or 2');
    }

    var out = [];
    C._uniqueTuples(f, [fwdCol, revCol]).forEach(function (pair) {
      var fwd = pair[0], rev = pair[1];
      var idx = [];
      f.rows.forEach(function (r, i) {
        if (r[fwdCol] === fwd && r[revCol] === rev) idx.push(i);
      });
      var first = f.rows[idx[0]];
      var tmCol = 'PCR' + n + ' Tm';
      if (!C.has(f, tmCol))
        throw new Error("The samples sheet has no '" + tmCol + "' column, " +
                        'which the thermocycler program is built from.');
      var tm = first[tmCol];
      var bpCol = C.has(f, 'PCR2 bp') ? 'PCR2 bp' : 'PCR1 bp';
      if (!C.has(f, bpCol))
        throw new Error("The samples sheet has no 'PCR1 bp' or 'PCR2 bp' " +
                        'column, which sets the extension time.');
      var bp = parseFloat(first[bpCol]);
      var ms = C._minSec(Math.floor(bp / 500) / 2 + 0.5);
      var annealTime;
      if (ms[0] === 0) annealTime = ms[1] + 's';
      else if (ms[1] === 0) annealTime = ms[0] + 'min';
      else annealTime = ms[0] + 'min ' + ms[1] + 's';

      var cyc = C.label({}, 'c', cycles) + ' cycles';
      out.push({key: fwd + '_' + rev + '_' + tm + '°C',
                table: {columns: ['Temperature', 'Time', 'Repeat'], rows: [
        {'Temperature': '98°C', 'Time': '30s', 'Repeat': ''},
        {'Temperature': '98°C', 'Time': '10s', 'Repeat': cyc},
        {'Temperature': tm + '°C', 'Time': '30s', 'Repeat': cyc},
        {'Temperature': '72°C', 'Time': annealTime, 'Repeat': cyc},
        {'Temperature': '72°C', 'Time': '2min', 'Repeat': ''},
        {'Temperature': '4°C', 'Time': '∞', 'Repeat': ''},
        {'Temperature': '', 'Time': '', 'Repeat': ''}]}});
    });
    return out;
  };

  /* thermocycler() assigns a whole column in place, which the Python does to
     the caller's frame. Mirrored, because the order the pcrs() body does
     things in depends on it. */
  C._setCol = function (f, col, value) {
    if (!C.has(f, col)) f.columns.push(col);
    f.rows.forEach(function (r) { r[col] = value; });
  };

  C._pivot = function (rows, values, indexCol, colCol) {
    var colVals = [];
    rows.forEach(function (r) {
      if (colVals.indexOf(r[colCol]) === -1) colVals.push(r[colCol]);
    });
    colVals.sort(function (a, b) { return a - b; });
    var keys = [], seen = {};
    rows.forEach(function (r) {
      var k = String(r[indexCol]) + '\u0000' + String(r['row']);
      if (!seen[k]) { seen[k] = true; keys.push([r[indexCol], r['row'], k]); }
    });
    keys.sort(function (x, y) {
      var a = x[0], b = y[0];
      if (a !== b) {
        if (typeof a === 'number' && typeof b === 'number') return a - b;
        return String(a) < String(b) ? -1 : 1;
      }
      return x[1] < y[1] ? -1 : (x[1] > y[1] ? 1 : 0);
    });
    var out = keys.map(function (k) {
      var o = {};
      colVals.forEach(function (c) {
        var hit = null;
        for (var i = 0; i < rows.length; i++) {
          var r = rows[i];
          if (r[indexCol] === k[0] && r['row'] === k[1] && r[colCol] === c) {
            hit = r[values];
            break;
          }
        }
        o[String(c)] = (hit === undefined) ? null : hit;
      });
      return o;
    });
    return {columns: colVals.map(String), rows: out};
  };

  C['edms.bio.ngs:pcrs'] = function (a) {
    var f = C.frame(a.df);
    // thermocycler() writes columns onto the frame it is given
    f = {columns: f.columns.slice(),
         rows: f.rows.map(function (r) { return C.merge(r, {}); })};

    var gdnaCol = a.gDNA_id_col || 'ID';
    var p1id = a.pcr1_id_col || 'PCR1 ID';
    var p1fwd = a.pcr1_fwd_col || 'PCR1 FWD';
    var p1rev = a.pcr1_rev_col || 'PCR1 REV';
    var p2id = a.pcr2_id_col || 'PCR2 ID';
    var p2fwd = a.pcr2_fwd_col || 'PCR2 FWD';
    var p2rev = a.pcr2_rev_col || 'PCR2 REV';
    var umiCol = a.umi_col || 'UMI';
    var pcr1TotalUL = C.g(a, 'pcr1_total_uL', 20);
    var pcr2TotalUL = C.g(a, 'pcr2_total_uL', 20);
    var mmX = C.g(a, 'mm_x', 1.1);
    var maxWellUL = C.g(a, 'max_well_uL', 50);
    var ultra = a.ultra === undefined ? false : Boolean(a.ultra);
    var syber = a.syber === undefined ? false : Boolean(a.syber);
    var inner1 = a.inner1 === undefined ? false : Boolean(a.inner1);
    var inner2 = a.inner2 === undefined ? false : Boolean(a.inner2);
    var pcr2 = a.pcr2 === undefined ? true : Boolean(a.pcr2);
    var split = a.split_pcr1_primers === undefined ? true
      : Boolean(a.split_pcr1_primers);
    var stocks = {
      Q5_mm_x_stock: C.g(a, 'Q5_mm_x_stock', 5),
      dNTP_mM_stock: C.g(a, 'dNTP_mM_stock', 10),
      fwd_uM_stock: C.g(a, 'fwd_uM_stock', 10),
      rev_uM_stock: C.g(a, 'rev_uM_stock', 10),
      Q5_U_uL_stock: C.g(a, 'Q5_U_uL_stock', 2),
      Q5_mm_x_desired: C.g(a, 'Q5_mm_x_desired', 1),
      dNTP_mM_desired: C.g(a, 'dNTP_mM_desired', 0.2),
      fwd_uM_desired: C.g(a, 'fwd_uM_desired', 0.5),
      rev_uM_desired: C.g(a, 'rev_uM_desired', 0.5),
      Q5_U_uL_desired: C.g(a, 'Q5_U_uL_desired', 0.02),
      syber_x_stock: 10, syber_x_desired: 1};

    // --- plate layouts ------------------------------------------------ //
    function plate(group, pcr) {
      var rowNames, colNums;
      if (group === '96-well plate' &&
          ((pcr === '1' && inner1) || (pcr === '2' && inner2))) {
        rowNames = ['B', 'C', 'D', 'E', 'F', 'G'];
        colNums = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
      } else if (group === '96-well plate') {
        rowNames = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'];
        colNums = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
      } else if (group === '8-strip plate') {
        rowNames = ['A', 'B', 'C', 'D', 'E', 'F'];
        colNums = [1, 2, 3, 4, 5, 6, 7, 8];
      } else {
        throw new Error("group must be '96-well plate' or '8-strip plate'");
      }
      var totalUL = pcr === '1' ? pcr1TotalUL : pcr2TotalUL;
      var wellULs = [];
      for (var k = 0; k < Math.trunc(totalUL / maxWellUL); k++)
        wellULs.push(maxWellUL);
      if (totalUL % maxWellUL > 0)
        wellULs.push(C.pyRound(totalUL % maxWellUL, 2));

      function fill(rows, plateI) {
        var out = [], rowI = 0, colI = 0;
        rows.forEach(function (src) {
          wellULs.forEach(function (uL) {
            if (colI >= colNums.length) {
              if (rowI >= rowNames.length - 1) { rowI = 0; colI = 0; plateI += 1; }
              else { rowI += 1; colI = 0; }
            }
            out.push(C.merge(src, {plate: plateI, row: rowNames[rowI],
                                   column: colNums[colI], uL: uL}));
            colI += 1;
          });
        });
        return {rows: out, plateI: plateI};
      }

      var placed;
      if (pcr === '1' && split) {
        placed = [];
        var plateI = 1;
        C._uniqueTuples(f, [p1fwd, p1rev]).forEach(function (pair) {
          var subset = f.rows.filter(function (r) {
            return r[p1fwd] === pair[0] && r[p1rev] === pair[1];
          });
          var res = fill(subset, plateI);
          res.rows.forEach(function (r) {
            r.plate = r.plate + '_' + pair[0] + '_' + pair[1];
          });
          placed = placed.concat(res.rows);
          plateI = res.plateI + 1;
        });
      } else {
        placed = fill(f.rows, 1).rows;
      }
      return placed;   // 'plate' is the index column, renamed by the caller
    }

    var wells96_1 = plate('96-well plate', '1');
    var wells96_2 = pcr2 ? plate('96-well plate', '2') : null;
    var strip1 = split ? plate('8-strip plate', '1') : null;
    var strip2 = (split && pcr2) ? plate('8-strip plate', '2') : null;

    var pivots = [];
    var add = function (name, table) { pivots.push({key: name, table: table}); };
    if (split) {
      add('96-well_' + gdnaCol, C._pivot(wells96_1, gdnaCol, 'plate', 'column'));
      add('96-well_' + p1id, C._pivot(wells96_1, p1id, 'plate', 'column'));
      if (pcr2)
        [p2id, p2fwd, p2rev].forEach(function (col) {
          add('96-well_' + col, C._pivot(wells96_2, col, 'plate', 'column'));
        });
      add('8-strip_' + gdnaCol, C._pivot(strip1, gdnaCol, 'plate', 'column'));
      add('8-strip_' + p1id, C._pivot(strip1, p1id, 'plate', 'column'));
      if (pcr2)
        [p2id, p2fwd, p2rev].forEach(function (col) {
          add('8-strip_' + col, C._pivot(strip2, col, 'plate', 'column'));
        });
    } else {
      [gdnaCol, p1id, p1fwd, p1rev].forEach(function (col) {
        add(col, C._pivot(wells96_1, col, 'plate', 'column'));
      });
      if (pcr2)
        [p2id, p2fwd, p2rev].forEach(function (col) {
          add(col, C._pivot(wells96_2, col, 'plate', 'column'));
        });
    }
    if (pcr1TotalUL > maxWellUL) {
      if (split) {
        add('96-well_PCR1 uL', C._pivot(wells96_1, 'uL', 'plate', 'column'));
        add('8-strip_PCR1 uL', C._pivot(strip1, 'uL', 'plate', 'column'));
      } else {
        add('PCR1 uL', C._pivot(wells96_1, 'uL', 'plate', 'column'));
      }
    }
    if (pcr2 && pcr2TotalUL > maxWellUL) {
      if (split) {
        add('96-well_PCR2 uL', C._pivot(wells96_2, 'uL', 'plate', 'column'));
        add('8-strip_PCR2 uL', C._pivot(strip2, 'uL', 'plate', 'column'));
      } else {
        add('PCR2 uL', C._pivot(wells96_2, 'uL', 'plate', 'column'));
      }
    }

    // --- master mixes and programs ------------------------------------ //
    var pcr1Primers = C._vcsOrdered(f, [p1fwd, p1rev]);
    var pcr2Primers = [];
    if (pcr2) {
      C._setCol(f, 'PCR2 FWD MM', 'PCR2-FWD');
      pcr2Primers = C._vcsOrdered(f, ['PCR2 FWD MM', p2rev]);
    }

    function mms(primers, template, templateUL, totalUL, withSyber) {
      var o = C.merge(stocks, {template: template, template_uL: templateUL,
                               total_uL: totalUL, mm_x: mmX,
                               syber: withSyber || false});
      if (!ultra) return C._pcrMM(primers, o);
      o.Q5_mm_x_stock = 2;   // NEBNext Ultra II Q5 master mix stock
      return C._pcrMMUltra(primers, o);
    }

    function ultraTemplateUL(totalUL, withSyber) {
      var parts = [stocks.Q5_mm_x_desired / 2,
                   stocks.fwd_uM_desired / stocks.fwd_uM_stock,
                   stocks.rev_uM_desired / stocks.rev_uM_stock];
      if (withSyber) parts.push(1 / 10);
      return totalUL - C._sumRep(parts, totalUL);
    }

    var hasUMI = C.has(f, umiCol);
    var pcr1Mms, pcr2Mms, thermos = [];
    if (!hasUMI) {
      pcr1Mms = mms(pcr1Primers, 'gDNA Extract',
                    ultra ? ultraTemplateUL(pcr1TotalUL, syber) : 5,
                    pcr1TotalUL, syber);
      pcr2Mms = pcr2 ? mms(pcr2Primers, 'PCR1 Product', 1, pcr2TotalUL, false)
                     : [];
      thermos = C._thermocycler(f, {n: '1', cycles: a.pcr1_cycles,
                                    pcr_fwd_col: p1fwd, pcr_rev_col: p1rev});
      if (pcr2)
        thermos = thermos.concat(
          C._thermocycler(f, {n: '2', cycles: a.pcr2_cycles,
                              pcr_fwd_col: p2fwd, pcr_rev_col: p2rev}));
    } else {
      var umiRows = f.rows.filter(function (r) { return r[umiCol] === true; });
      var nonUmiRows = f.rows.filter(function (r) { return r[umiCol] === false; });
      var fUmi = {columns: f.columns.slice(),
                  rows: umiRows.map(function (r) { return C.merge(r, {}); })};
      var fNon = {columns: f.columns.slice(),
                  rows: nonUmiRows.map(function (r) { return C.merge(r, {}); })};
      var tm15 = (a.pcr1_5_Tm === null || a.pcr1_5_Tm === undefined ||
                  a.pcr1_5_Tm === '') ? 65 : a.pcr1_5_Tm;
      fUmi.rows.forEach(function (r) { r['PCR1.5 ID'] = r[p1id] + '_UMI'; });
      ['PCR1.5 ID', 'PCR1.5 FWD', 'PCR1.5 REV', 'PCR1.5 Tm']
        .forEach(function (c) {
          if (fUmi.columns.indexOf(c) === -1) fUmi.columns.push(c);
        });
      C._setCol(fUmi, 'PCR1.5 FWD', 'P5-FWD');
      C._setCol(fUmi, 'PCR1.5 REV', 'P7-REV');
      C._setCol(fUmi, 'PCR1.5 Tm', tm15);

      pcr1Mms = mms(pcr1Primers, 'gDNA Extract',
                    ultra ? ultraTemplateUL(pcr1TotalUL, syber) : 5,
                    pcr1TotalUL, syber);
      /* edms.bio.ngs:pcrs builds a PCR1.5 master mix here and then leaves it
         out of what it returns, so it never reaches a rendered protocol.
         Not built here either, for the same reason: this file's job is to
         agree with the package. */
      pcr2Mms = pcr2 ? mms(pcr2Primers, 'PCR1.5 Product', 1, pcr2TotalUL, false)
                     : [];
      thermos = C._thermocycler(fUmi, {n: '1', cycles: a.umi_cycles,
                                       pcr_fwd_col: p1fwd, pcr_rev_col: p1rev,
                                       umi: true})
        .concat(C._thermocycler(fUmi, {n: '1.5', cycles: a.pcr1_cycles,
                                       pcr_fwd_col: 'PCR1.5 FWD',
                                       pcr_rev_col: 'PCR1.5 REV', umi: true}))
        .concat(C._thermocycler(fNon, {n: '1', cycles: a.pcr1_cycles,
                                       pcr_fwd_col: p1fwd,
                                       pcr_rev_col: p1rev}));
      if (pcr2)
        thermos = thermos.concat(
          C._thermocycler(f, {n: '2', cycles: a.pcr2_cycles,
                              pcr_fwd_col: p2fwd, pcr_rev_col: p2rev}));
    }

    /* edms.work.plan._tables walks the returned tuple of dicts and names each
       table after its key, disambiguating a repeat with _2. */
    var tables = {};
    var put = function (key, table) {
      var name = String(key);
      if (Object.prototype.hasOwnProperty.call(tables, name)) {
        var i = 2;
        while (Object.prototype.hasOwnProperty.call(tables, name + '_' + i)) i++;
        name = name + '_' + i;
      }
      tables[name] = table;
    };
    pivots.forEach(function (p) { put(p.key, p.table); });
    pcr1Mms.forEach(function (m) { put(C.pyTuple(m.key), m.table); });
    pcr2Mms.forEach(function (m) { put(C.pyTuple(m.key), m.table); });
    thermos.forEach(function (t) { put(t.key, t.table); });
    return {__tables__: tables};
  };

})(EDMS.calc);

if (typeof module !== 'undefined' && module.exports) module.exports = EDMS;
