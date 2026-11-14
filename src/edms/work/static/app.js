/*
 * app.js -- make a rendered protocol recalculate in place.
 *
 * Editing a parameter here changes nothing on disk. The page recomputes the
 * tables whose calculation calc.js implements, marks itself unsaved, and can
 * hand back the edms work deviate command that would make the change part of
 * the record. A protocol is a record; a scratch calculation that silently
 * became one would be worse than no calculator at all.
 */
'use strict';

(function () {
  var el = document.getElementById('edms-data');
  if (!el || !window.EDMS || !EDMS.calc) return;

  var D;
  try { D = JSON.parse(el.textContent); } catch (e) { return; }
  EDMS.consts = D.consts;

  var base = JSON.parse(JSON.stringify(D.params));
  var params = JSON.parse(JSON.stringify(D.params));
  var inputs = JSON.parse(JSON.stringify(D.inputs || {}));

  // ------------------------------------------------------------------ //
  // Argument resolution, mirroring plan._resolve_arg
  // ------------------------------------------------------------------ //
  function resolveArgs(spec) {
    var out = {};
    Object.keys(spec || {}).forEach(function (k) {
      var v = spec[k];
      if (typeof v === 'string' && v.charAt(0) === '@') out[k] = inputs[v.slice(1)];
      else if (typeof v === 'string' && v.charAt(0) === '$') out[k] = params[v.slice(1)];
      else out[k] = v;
    });
    return out;
  }

  /* Which steps use a given parameter: that is the step a deviation on it
     belongs to, the same way edms work deviate takes one. */
  function stepsUsing(name) {
    return D.steps.filter(function (s) {
      return Object.keys(s.args || {}).some(function (k) {
        return s.args[k] === '$' + name;
      });
    }).map(function (s) { return s.id; });
  }

  // ------------------------------------------------------------------ //
  // Tables
  // ------------------------------------------------------------------ //
  function fmtCell(v) {
    if (v === null || v === undefined) return '';
    if (typeof v === 'number') {
      if (Number.isInteger(v)) return String(v);
      return String(Math.abs(v) < 1e-4 ? v.toExponential(2) : v);
    }
    return String(v);
  }

  function paint(host, table) {
    var h = ['<table><thead><tr>'];
    table.columns.forEach(function (c) { h.push('<th>' + EDMS.plate.esc(c) + '</th>'); });
    h.push('</tr></thead><tbody>');
    if (!table.rows.length)
      h.push('<tr><td class="pv-empty" colspan="' + table.columns.length +
             '">nothing here &mdash; which is the good case</td></tr>');
    table.rows.forEach(function (r) {
      h.push('<tr>');
      table.columns.forEach(function (c) {
        h.push('<td>' + EDMS.plate.esc(fmtCell(r[c])) + '</td>');
      });
      h.push('</tr>');
    });
    h.push('</tbody></table>');
    host.innerHTML = h.join('');
  }

  // ------------------------------------------------------------------ //
  // Plate map
  // ------------------------------------------------------------------ //
  function conditionLabels(nCond) {
    var plasmids = inputs.plasmids;
    if (plasmids && plasmids.rows) {
      var vs = params.VSVG_plasmid, gp = params.GagPol_plasmid;
      var transfers = plasmids.rows.filter(function (r) {
        return r.Plasmid !== vs && r.Plasmid !== gp;
      }).map(function (r) { return String(r.Plasmid); });
      if (transfers.length) return transfers;
    }
    var out = [];
    for (var i = 0; i < nCond; i++) out.push('Condition ' + (i + 1));
    return out;
  }

  function drawPlate(host, geomRow, fmt) {
    var used = geomRow && geomRow['Wells needed'];
    var reps = (geomRow && geomRow['Replicates']) || 1;
    var conds = (geomRow && geomRow['Conditions']) || 1;
    var spec = EDMS.calc.spec(fmt);
    if (!used) used = spec.wells === 1 ? 1 : conds * reps;

    var names = conditionLabels(conds);
    var labels = [], groups = [];
    for (var c = 0; c < conds; c++) {
      for (var r = 0; r < reps; r++) {
        var nm = names[c] || ('Condition ' + (c + 1));
        labels.push(nm + (reps > 1 ? ' rep ' + (r + 1) : ''));
        groups.push(nm);
      }
    }
    var perPlate = spec.wells;
    var plates = Math.max(1, Math.ceil(used / perPlate));
    var svgs = [];
    for (var p = 0; p < plates; p++) {
      var from = p * perPlate;
      svgs.push(EDMS.plate.svg({
        format: fmt,
        used: Math.max(0, Math.min(used - from, perPlate)),
        labels: labels.slice(from),
        groups: groups.slice(from)
      }));
    }
    host.innerHTML = '<div class="pv-wrap">' + svgs.join('') + '</div>' +
      EDMS.plate.legend(groups);
  }

  // ------------------------------------------------------------------ //
  // Recompute
  // ------------------------------------------------------------------ //
  var errorsBox = null;

  function recompute() {
    var errs = [];
    D.steps.forEach(function (s) {
      if (!s.live) return;
      var hosts = document.querySelectorAll('[data-step="' + s.id + '"]');
      if (!hosts.length) return;
      var result;
      try {
        result = EDMS.calc[s.compute](resolveArgs(s.args));
      } catch (e) {
        errs.push(s.id + ': ' + e.message);
        Array.prototype.forEach.call(hosts, function (h) {
          h.innerHTML = '<div class="err">Could not compute: ' +
            EDMS.plate.esc(e.message) + '</div>';
        });
        return;
      }
      var tables = result.__tables__ || null;
      Array.prototype.forEach.call(hosts, function (h, i) {
        paint(h, tables ? tables[h.getAttribute('data-table')] : result);
      });
      if (s.compute === 'edms.work.plates:geometry') {
        var map = document.querySelector('[data-plate="' + s.id + '"]');
        if (map) drawPlate(map, result.rows[0], resolveArgs(s.args).plate_format);
      }
    });
    if (errorsBox)
      errorsBox.innerHTML = errs.length
        ? '<div class="err">' + errs.map(EDMS.plate.esc).join('<br>') + '</div>' : '';
    markDirty();
  }

  // ------------------------------------------------------------------ //
  // Controls
  // ------------------------------------------------------------------ //
  function changed() {
    return Object.keys(params).filter(function (k) {
      return String(params[k]) !== String(base[k]);
    });
  }

  function markDirty() {
    var diff = changed();
    var bar = document.getElementById('edms-dirty');
    if (!bar) return;
    if (!diff.length) {
      bar.className = 'edms-dirty';
      bar.innerHTML = 'Showing the run as recorded.';
      return;
    }
    bar.className = 'edms-dirty on';
    bar.innerHTML = '<b>Scratch:</b> ' + diff.map(function (k) {
      return EDMS.plate.esc(k) + ' ' + EDMS.plate.esc(base[k]) + ' &rarr; ' +
        EDMS.plate.esc(params[k]);
    }).join(', ') + '. Nothing on disk has changed.' +
      '<button type="button" id="edms-copy">Copy deviate command' +
      (diff.length > 1 ? 's' : '') + '</button>' +
      '<button type="button" id="edms-reset">Reset</button>';
    document.getElementById('edms-reset').onclick = function () {
      params = JSON.parse(JSON.stringify(base));
      buildPanel();
      recompute();
    };
    document.getElementById('edms-copy').onclick = function () {
      var why = (document.getElementById('edms-why') || {}).value || '';
      var cmds = diff.map(function (k) {
        var step = stepsUsing(k)[0] || (D.steps[0] || {}).id || '<step>';
        return 'edms work deviate ' + D.run_id + ' -s ' + step + ' -p ' + k +
          ' -v "' + params[k] + '"' + (why ? ' -w "' + why + '"' : '');
      }).join('\n');
      navigator.clipboard.writeText(cmds).then(function () {
        var b = document.getElementById('edms-copy');
        b.textContent = 'Copied — paste in the terminal';
      }, function () {
        window.prompt('Copy this:', cmds);
      });
    };
  }

  function control(name, value) {
    if (name === 'plate_format') {
      var sel = document.createElement('select');
      D.formats.forEach(function (f) {
        var o = document.createElement('option');
        o.value = f; o.textContent = f;
        if (f === value) o.selected = true;
        sel.appendChild(o);
      });
      return sel;
    }
    var inp = document.createElement('input');
    if (typeof value === 'number') {
      inp.type = 'number';
      inp.step = Number.isInteger(value) ? '1' : 'any';
    } else {
      inp.type = 'text';
    }
    inp.value = value;
    return inp;
  }

  function buildPanel() {
    var host = document.getElementById('edms-controls');
    if (!host) return;
    host.innerHTML = '';
    Object.keys(params).forEach(function (k) {
      var used = stepsUsing(k);
      var wrap = document.createElement('label');
      wrap.className = 'edms-ctl' + (used.length ? '' : ' inert');
      var lab = document.createElement('span');
      lab.textContent = k;
      if (!used.length) lab.title = 'Used in the step text only, not in a calculation';
      var c = control(k, params[k]);
      c.addEventListener('change', function () {
        var v = c.value;
        params[k] = (typeof base[k] === 'number' && v !== '' && isFinite(+v))
          ? +v : v;
        recompute();
      });
      wrap.appendChild(lab);
      wrap.appendChild(c);
      host.appendChild(wrap);
    });
  }

  function makeInputsEditable() {
    Array.prototype.forEach.call(
      document.querySelectorAll('[data-input]'), function (host) {
        var name = host.getAttribute('data-input');
        var data = inputs[name];
        if (!data) return;
        var tbl = host.querySelector('table');
        if (!tbl) return;
        var bodyRows = tbl.querySelectorAll('tbody tr');
        Array.prototype.forEach.call(bodyRows, function (tr, ri) {
          Array.prototype.forEach.call(tr.children, function (td, ci) {
            var col = data.columns[ci];
            if (typeof data.rows[ri][col] !== 'number') return;
            td.contentEditable = 'true';
            td.className = 'edit';
            td.addEventListener('blur', function () {
              var v = parseFloat(td.textContent);
              if (isFinite(v)) {
                data.rows[ri][col] = v;
                recompute();
              } else {
                td.textContent = data.rows[ri][col];
              }
            });
            td.addEventListener('keydown', function (ev) {
              if (ev.key === 'Enter') { ev.preventDefault(); td.blur(); }
            });
          });
        });
      });
  }

  // ------------------------------------------------------------------ //
  function start() {
    errorsBox = document.getElementById('edms-errors');
    buildPanel();
    makeInputsEditable();
    recompute();
    var notLive = D.steps.filter(function (s) { return !s.live; });
    var warn = document.getElementById('edms-static');
    if (warn && notLive.length)
      warn.innerHTML = 'Static here: ' + notLive.map(function (s) {
        return '<code>' + EDMS.plate.esc(s.id) + '</code>';
      }).join(', ') + '. Those recalculate with <code>edms work plan</code>.';
  }

  if (document.readyState === 'loading')
    document.addEventListener('DOMContentLoaded', start);
  else start();
})();
