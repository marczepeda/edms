/*
 * plate.js -- draw the vessel a protocol is run in.
 *
 * A protocol says '6 wells of a 96-well plate' and leaves you to picture it.
 * This draws it, so the shape of the experiment is visible before you are
 * standing over the plate with a multichannel.
 *
 * Grid shapes are keyed by well count, which is layout and not arithmetic:
 * the areas and media volumes still come from the Python in EDMS.consts.
 */
'use strict';

var EDMS = (typeof EDMS !== 'undefined') ? EDMS : {};
EDMS.plate = {};

(function (P) {

  var GRID = {96: [8, 12], 48: [6, 8], 24: [4, 6], 12: [3, 4], 6: [2, 3]};
  var ROWS = 'ABCDEFGH';

  /* Distinguishable at small sizes and in both themes. Not a gradient: these
     label conditions, which have no order. */
  P.PALETTE = ['#4f8ff7', '#e8833a', '#3fa66a', '#c3539b',
               '#d9b43c', '#6f6bd8', '#4bb3c4', '#cc5f5f'];

  P.wellName = function (i, cols) {
    return ROWS[Math.floor(i / cols)] + String((i % cols) + 1);
  };

  /*
   * svg(): the vessel, with the used wells filled.
   *
   * opts.format   canonical format name
   * opts.used     how many wells are in use
   * opts.wells    the well names in use ('B2', 'C7', ...), when they are not
   *               simply the first opts.used of them -- a plate filled with
   *               the edge wells left out is not filled from A1
   * opts.labels   one label per used well, aligned to wells or fill order
   * opts.groups   one group key per used well, colours wells together
   */
  P.svg = function (opts) {
    var spec = EDMS.calc.spec(opts.format);
    var n = spec.wells;
    var used = Math.max(0, Math.min(opts.used === undefined ? 0 : opts.used, n));
    var labels = opts.labels || [];
    var groups = opts.groups || [];

    var keys = [];
    groups.forEach(function (g) { if (keys.indexOf(g) === -1) keys.push(g); });
    var colourOf = function (i) {
      if (!groups.length) return P.PALETTE[0];
      return P.PALETTE[keys.indexOf(groups[i]) % P.PALETTE.length];
    };

    var parts, w, h;
    if (!GRID[n]) {
      // A dish: one vessel, drawn to scale against the largest we know
      var r = 46;
      w = 150; h = 150;
      parts = ['<circle cx="75" cy="72" r="' + r + '" class="pv-dish' +
               (used > 0 ? ' pv-on' : '') + '"' +
               (used > 0 ? ' style="fill:' + colourOf(0) + '"' : '') + '/>'];
      if (labels[0])
        parts.push('<text x="75" y="76" class="pv-lab pv-lab-dish">' +
                   P.esc(labels[0]) + '</text>');
      parts.push('<text x="75" y="138" class="pv-cap">' +
                 P.esc(spec.plate_format) + ' &middot; ' + spec.area_cm2 +
                 ' cm&sup2;</text>');
    } else {
      var rows = GRID[n][0], cols = GRID[n][1];
      /* Which well holds which sample. Without opts.wells this is the first
         opts.used wells in fill order, which is what a plan that fills from
         A1 means; with it, a well is on only if it was named. */
      var slotOf = {};
      if (opts.wells && opts.wells.length) {
        opts.wells.forEach(function (name, k) {
          var r = ROWS.indexOf(String(name).charAt(0));
          var c = parseInt(String(name).slice(1), 10) - 1;
          if (r >= 0 && c >= 0 && r < rows && c < cols) slotOf[r * cols + c] = k;
        });
      } else {
        for (var u = 0; u < used; u++) slotOf[u] = u;
      }
      var pitch = cols > 8 ? 26 : 34;
      var rad = pitch * 0.40;
      var padL = 22, padT = 20;
      w = padL + cols * pitch + 10;
      h = padT + rows * pitch + 26;
      parts = ['<rect x="' + (padL - 10) + '" y="' + (padT - 10) + '" rx="6" ' +
               'width="' + (cols * pitch + 8) + '" height="' + (rows * pitch + 8) +
               '" class="pv-body"/>'];
      for (var c = 0; c < cols; c++)
        parts.push('<text x="' + (padL + c * pitch + pitch / 2 - 1) + '" y="' +
                   (padT - 1) + '" class="pv-ax">' + (c + 1) + '</text>');
      for (var r2 = 0; r2 < rows; r2++)
        parts.push('<text x="' + (padL - 14) + '" y="' +
                   (padT + r2 * pitch + pitch / 2 + 4) + '" class="pv-ax">' +
                   ROWS[r2] + '</text>');
      var filled = 0;
      for (var i = 0; i < n; i++) {
        var cx = padL + (i % cols) * pitch + pitch / 2;
        var cy = padT + Math.floor(i / cols) * pitch + pitch / 2;
        var k = slotOf[i];
        var on = k !== undefined;
        if (on) filled += 1;
        parts.push('<circle cx="' + cx + '" cy="' + cy + '" r="' + rad +
                   '" class="pv-well' + (on ? ' pv-on' : '') + '"' +
                   (on ? ' style="fill:' + colourOf(k) + '"' : '') +
                   '><title>' + P.esc(P.wellName(i, cols) +
                   (on && labels[k] ? ' - ' + labels[k] : '')) +
                   '</title></circle>');
      }
      parts.push('<text x="' + (padL - 10) + '" y="' + (h - 6) +
                 '" class="pv-cap pv-cap-left">' + P.esc(spec.plate_format) +
                 ' &middot; ' + filled + ' of ' + n + ' wells &middot; ' +
                 spec.media_mL + ' mL per well</text>');
      used = filled;
    }
    return '<svg class="pv" viewBox="0 0 ' + w + ' ' + h + '" width="' + w +
           '" height="' + h + '" role="img" aria-label="' +
           P.esc(spec.plate_format + ', ' + used + ' of ' + n + ' wells in use') +
           '">' + parts.join('') + '</svg>';
  };

  /* legend(): a swatch per group key, in the order the wells use them */
  P.legend = function (groups) {
    var keys = [];
    (groups || []).forEach(function (g) { if (keys.indexOf(g) === -1) keys.push(g); });
    if (keys.length < 2) return '';
    return '<div class="pv-key">' + keys.map(function (k, i) {
      return '<span><i style="background:' + P.PALETTE[i % P.PALETTE.length] +
             '"></i>' + P.esc(k) + '</span>';
    }).join('') + '</div>';
  };

  P.esc = function (s) {
    return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;')
      .replace(/>/g, '&gt;').replace(/"/g, '&quot;');
  };

})(EDMS.plate);

if (typeof module !== 'undefined' && module.exports) module.exports = EDMS;
