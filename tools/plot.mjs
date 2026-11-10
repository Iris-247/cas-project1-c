// Optional figure generator for the C project.  Reads the .dat tables in
// output/ and writes PNGs to figures/.  Requires Node >= 18 and
// @napi-rs/canvas (npm install in this directory).
//
//   node tools/plot.mjs            # output/ -> figures/
//
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { createCanvas } from "@napi-rs/canvas";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const DATA = path.join(ROOT, "output");
const FIGS = path.join(ROOT, "figures");
fs.mkdirSync(FIGS, { recursive: true });

const W = 1400, H = 900;
const L = 130, R = W - 55, T = 85, B = H - 105;
const FONT = '"Segoe UI", Arial, sans-serif';

function readDat(name) {
    const file = path.join(DATA, name);
    const lines = fs.readFileSync(file, "utf8").split(/\r?\n/);
    const rows = [];
    for (const line of lines) {
        const s = line.trim();
        if (!s || s.startsWith("#")) continue;
        const vals = s.split(/\s+/).map(Number);
        if (vals.every(Number.isFinite)) rows.push(vals);
    }
    return rows;
}

function logTicks(min, max) {
    const out = [];
    for (let e = Math.ceil(Math.log10(min)); e <= Math.floor(Math.log10(max)); e++)
        out.push(Math.pow(10, e));
    return out;
}

function linTicks(min, max) {
    const span = max - min;
    const raw = span / 5;
    const p = Math.pow(10, Math.floor(Math.log10(raw)));
    let step = p;
    for (const m of [1, 2, 2.5, 5, 10]) {
        if (p * m >= raw) { step = p * m; break; }
    }
    const out = [];
    for (let v = Math.ceil(min / step) * step; v <= max + step * 1e-9; v += step)
        out.push(v);
    return out;
}

function fmt(v, log) {
    if (log) {
        const e = Math.round(Math.log10(v));
        return `1e${e}`;
    }
    if (v === 0) return "0";
    const a = Math.abs(v);
    if (a >= 1e4 || a < 1e-2) return v.toExponential(1).replace("e+", "e");
    return Number(v.toPrecision(3)).toString();
}

function drawFigure(o) {
    const canvas = createCanvas(W, H);
    const ctx = canvas.getContext("2d");
    const [xmin, xmax] = o.xlim;
    const [ymin, ymax] = o.ylim;
    const xlog = !!o.xlog, ylog = !!o.ylog;

    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, W, H);

    const mx = (v) => {
        const t = xlog
            ? (Math.log10(v) - Math.log10(xmin)) / (Math.log10(xmax) - Math.log10(xmin))
            : (v - xmin) / (xmax - xmin);
        return L + t * (R - L);
    };
    const my = (v) => {
        const t = ylog
            ? (Math.log10(v) - Math.log10(ymin)) / (Math.log10(ymax) - Math.log10(ymin))
            : (v - ymin) / (ymax - ymin);
        return B - t * (B - T);
    };

    ctx.strokeStyle = "#111";
    ctx.lineWidth = 2.2;
    ctx.strokeRect(L, T, R - L, B - T);

    // ---- ticks ----
    ctx.font = `19px ${FONT}`;
    ctx.fillStyle = "#111";
    const xt = xlog ? logTicks(xmin, xmax) : linTicks(xmin, xmax);
    const yt = ylog ? logTicks(ymin, ymax) : linTicks(ymin, ymax);

    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    for (const v of xt) {
        const x = mx(v);
        ctx.beginPath(); ctx.moveTo(x, B); ctx.lineTo(x, B + 8); ctx.stroke();
        ctx.fillText(fmt(v, xlog), x, B + 14);
        if (xlog)
            for (let m = 2; m <= 9; m++) {
                const vm = v * m;
                if (vm >= xmax) break;
                const xm = mx(vm);
                ctx.beginPath(); ctx.moveTo(xm, B); ctx.lineTo(xm, B + 4); ctx.stroke();
            }
    }

    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    for (const v of yt) {
        const y = my(v);
        ctx.beginPath(); ctx.moveTo(L, y); ctx.lineTo(L - 8, y); ctx.stroke();
        ctx.fillText(fmt(v, ylog), L - 14, y);
        if (ylog)
            for (let m = 2; m <= 9; m++) {
                const vm = v * m;
                if (vm >= ymax) break;
                const ym = my(vm);
                ctx.beginPath(); ctx.moveTo(L, ym); ctx.lineTo(L - 4, ym); ctx.stroke();
            }
    }

    // ---- labels / title ----
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    ctx.font = `22px ${FONT}`;
    ctx.fillText(o.xlabel, (L + R) / 2, H - 55);
    ctx.save();
    ctx.translate(42, (T + B) / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textBaseline = "middle";
    ctx.fillText(o.ylabel, 0, 0);
    ctx.restore();

    ctx.font = `bold 26px ${FONT}`;
    ctx.textBaseline = "top";
    ctx.fillText(o.title, (L + R) / 2, 28);

    // ---- series ----
    ctx.save();
    ctx.beginPath();
    ctx.rect(L, T, R - L, B - T);
    ctx.clip();

    for (const s of o.series) {
        if (s.hline !== undefined) {
            const y = my(s.hline);
            ctx.setLineDash(s.dash || []);
            ctx.strokeStyle = s.color;
            ctx.lineWidth = s.width || 2.2;
            ctx.beginPath(); ctx.moveTo(L, y); ctx.lineTo(R, y); ctx.stroke();
            continue;
        }
        ctx.setLineDash(s.dash || []);
        ctx.strokeStyle = s.color;
        ctx.lineWidth = s.width || 2.4;
        ctx.beginPath();
        let started = false;
        for (const row of s.data) {
            const xv = row[s.xcol ?? 0], yv = row[s.ycol ?? 1];
            if (!Number.isFinite(xv) || !Number.isFinite(yv)) continue;
            if (xlog && xv <= 0) continue;
            if (ylog && yv <= 0) continue;
            if (xv < xmin || xv > xmax || yv < ymin || yv > ymax) {
                started = false;
                continue;
            }
            const px = mx(xv), py = my(yv);
            if (!started) { ctx.moveTo(px, py); started = true; }
            else ctx.lineTo(px, py);
        }
        ctx.stroke();
    }
    ctx.restore();
    ctx.setLineDash([]);

    // ---- legend ----
    if (o.series.some((s) => s.label)) {
        ctx.font = `19px ${FONT}`;
        let w = 0;
        for (const s of o.series) if (s.label) w = Math.max(w, ctx.measureText(s.label).width);
        const boxW = w + 90;
        const boxH = 34 * o.series.filter((s) => s.label).length + 16;
        const pos = o.legend || "tr";
        const bx = (pos === "tl" || pos === "bl") ? L + 18 : R - boxW - 18;
        const by = (pos === "tr" || pos === "tl") ? T + 18 : B - boxH - 18;

        ctx.fillStyle = "rgba(255,255,255,0.92)";
        ctx.strokeStyle = "#999";
        ctx.lineWidth = 1.2;
        ctx.fillRect(bx, by, boxW, boxH);
        ctx.strokeRect(bx, by, boxW, boxH);

        let yy = by + 16;
        for (const s of o.series) {
            if (!s.label) continue;
            ctx.strokeStyle = s.color;
            ctx.lineWidth = s.width || 2.4;
            ctx.setLineDash(s.dash || []);
            ctx.beginPath(); ctx.moveTo(bx + 14, yy + 9); ctx.lineTo(bx + 54, yy + 9); ctx.stroke();
            ctx.setLineDash([]);
            ctx.fillStyle = "#111";
            ctx.textAlign = "left";
            ctx.textBaseline = "middle";
            ctx.fillText(s.label, bx + 64, yy + 10);
            yy += 34;
        }
    }

    return canvas;
}

function save(canvas, name) {
    const file = path.join(FIGS, name);
    fs.writeFileSync(file, canvas.toBuffer("image/png"));
    console.log("wrote", path.relative(ROOT, file));
}

function filterX(rows, xmin) {
    return rows.filter((r) => r[0] >= xmin);
}

function main() {
    // 1. dark-matter mass (+ BCG for context)
    {
        const m = readDat("mass.dat");
        // Hernquist BCG mass, computed here so mass.dat stays a DM-only table
        const mbcg = m.map((r) => {
            const a = 4.971; // kpc
            return [r[0], 1.0e12 * r[0] * r[0] / ((r[0] + a) * (r[0] + a))];
        });
        save(drawFigure({
            title: "Dark matter and BCG mass profiles",
            xlabel: "r [kpc]", ylabel: "M(<r) [M_sun]",
            xlim: [0.1, 3000], ylim: [1e7, 1e16], xlog: true, ylog: true,
            series: [
                { data: m, xcol: 0, ycol: 1, color: "#111", label: "DM shell sum" },
                { data: m, xcol: 0, ycol: 2, color: "#d62728", dash: [10, 6], label: "NFW analytic" },
                { data: mbcg, xcol: 0, ycol: 1, color: "#2ca02c", dash: [4, 4], label: "BCG (Hernquist)" },
            ],
            legend: "tl",
        }), "mass.png");
    }

    // 2. isothermal density vs analytic
    {
        const g = readDat("gas.dat");
        save(drawFigure({
            title: "Isothermal gas density (DM only)",
            xlabel: "r [kpc]", ylabel: "rho [g cm^-3]",
            xlim: [0.1, 3000], ylim: [1e-30, 1e-22], xlog: true, ylog: true,
            series: [
                { data: g, xcol: 0, ycol: 1, color: "#111", label: "numerical HSE" },
                { data: g, xcol: 0, ycol: 2, color: "#d62728", dash: [10, 6], label: "analytic" },
                { data: g, xcol: 0, ycol: 3, color: "#1f77b4", dash: [3, 4], label: "NFW" },
            ],
            legend: "tr",
        }), "density_iso.png");
    }

    // 3. all density models
    {
        const iso = readDat("gas.dat");
        const bcg = readDat("gas_BCG.dat");
        const non = readDat("gas_BCG_non-iso.dat");
        const tur = readDat("gas_BCG_non-iso-turb.dat");
        save(drawFigure({
            title: "Gas density models",
            xlabel: "r [kpc]", ylabel: "rho [g cm^-3]",
            xlim: [0.1, 3000], ylim: [1e-30, 1e-22], xlog: true, ylog: true,
            series: [
                { data: iso, xcol: 0, ycol: 1, color: "#1f77b4", label: "isothermal, no BCG" },
                { data: bcg, xcol: 0, ycol: 1, color: "#2ca02c", label: "isothermal + BCG" },
                { data: non, xcol: 0, ycol: 1, color: "#111", label: "T(r) + BCG" },
                { data: tur, xcol: 0, ycol: 1, color: "#d62728", label: "T(r) + BCG + P_nt" },
                { data: iso, xcol: 0, ycol: 3, color: "#999", dash: [2, 5], label: "NFW" },
            ],
            legend: "tr",
        }), "density_models.png");
    }

    // 4. temperature
    {
        const t = readDat("temperature.dat");
        save(drawFigure({
            title: "Vikhlinin-type temperature profile",
            xlabel: "r [kpc]", ylabel: "T [K]",
            xlim: [0.1, 3000], ylim: [3e7, 1.2e8], xlog: true, ylog: false,
            series: [
                { data: t, xcol: 0, ycol: 1, color: "#111", label: "T(r)" },
                { hline: 8.9e7, color: "#d62728", dash: [10, 6], label: "T_mg" },
            ],
            legend: "tl",
        }), "temperature.png");
    }

    // 5. turbulence
    {
        const tur = readDat("gas_BCG_non-iso-turb.dat");
        save(drawFigure({
            title: "Effect of non-thermal pressure (v_t = 500 km/s)",
            xlabel: "r [kpc]", ylabel: "rho [g cm^-3]",
            xlim: [0.1, 3000], ylim: [1e-30, 1e-22], xlog: true, ylog: true,
            series: [
                { data: tur, xcol: 0, ycol: 2, color: "#111", label: "T(r), no P_nt" },
                { data: tur, xcol: 0, ycol: 1, color: "#d62728", label: "T(r) + P_nt" },
            ],
            legend: "tr",
        }), "turbulence.png");
    }

    // 6. Fe diffusion only
    {
        const f = readDat("fe_diffusion.dat");
        save(drawFigure({
            title: "Fe diffusion only (initial observed profile - Z_out)",
            xlabel: "r [kpc]", ylabel: "Z_Fe [Z_sun]",
            xlim: [0.1, 3000], ylim: [0, 0.6], xlog: true, ylog: false,
            series: [
                { data: f, xcol: 0, ycol: 1, color: "#111", label: "t = 0" },
                { data: f, xcol: 0, ycol: 2, color: "#1f77b4", label: "t = 1 Gyr" },
                { data: f, xcol: 0, ycol: 3, color: "#2ca02c", label: "t = 2 Gyr" },
                { data: f, xcol: 0, ycol: 4, color: "#d62728", label: "t = 5 Gyr" },
            ],
        }), "fe_diffusion.png");
    }

    // 7. Fe source only (clip the unphysically large core values)
    {
        const f = filterX(readDat("fe_source.dat"), 1.0);
        save(drawFigure({
            title: "Fe source only (Z_Fe = 0 initially)",
            xlabel: "r [kpc]", ylabel: "Z_Fe [Z_sun]",
            xlim: [1, 3000], ylim: [0, 5], xlog: true, ylog: false,
            series: [
                { data: f, xcol: 0, ycol: 1, color: "#1f77b4", label: "t = 1 Gyr" },
                { data: f, xcol: 0, ycol: 2, color: "#d62728", label: "t = 5 Gyr" },
            ],
        }), "fe_source.png");
    }

    // 8. source + diffusion vs observed excess
    {
        const f = readDat("fe_both.dat");
        const obs = readDat("fe_observed.dat");
        save(drawFigure({
            title: "Fe source + diffusion vs observed excess",
            xlabel: "r [kpc]", ylabel: "Z_Fe [Z_sun]",
            xlim: [0.1, 3000], ylim: [0, 0.7], xlog: true, ylog: false,
            series: [
                { data: obs, xcol: 0, ycol: 2, color: "#111", dash: [4, 4], label: "observed - Z_out" },
                { data: f, xcol: 0, ycol: 1, color: "#1f77b4", label: "model, t = 1 Gyr" },
                { data: f, xcol: 0, ycol: 2, color: "#d62728", label: "model, t = 5 Gyr" },
            ],
        }), "fe_both.png");
    }

    // 9. Fe mass budget
    {
        // fe_masses.dat has a text column, so parse it by hand.
        const txt = fs.readFileSync(path.join(DATA, "fe_masses.dat"), "utf8");
        const sets = { "diffusion only": [], "source only": [], "source + diffusion": [] };
        for (const line of txt.split(/\r?\n/)) {
            const s = line.trim();
            if (!s || s.startsWith("#")) continue;
            const p = s.split(/\s+/);
            const t = Number(p[0]), name = p[1], m100 = Number(p[2]);
            if (!Number.isFinite(t) || !Number.isFinite(m100)) continue;
            const key = name === "diffusion_only" ? "diffusion only"
                      : name === "source_only" ? "source only" : "source + diffusion";
            sets[key].push([t, m100]);
        }
        save(drawFigure({
            title: "Iron mass within 100 kpc",
            xlabel: "t [Gyr]", ylabel: "M_Fe(<100 kpc) [M_sun]",
            xlim: [0, 5], ylim: [1e7, 1e9], xlog: false, ylog: true,
            series: [
                { data: sets["diffusion only"], xcol: 0, ycol: 1, color: "#1f77b4", label: "diffusion only" },
                { data: sets["source only"], xcol: 0, ycol: 1, color: "#2ca02c", label: "source only" },
                { data: sets["source + diffusion"], xcol: 0, ycol: 1, color: "#d62728", label: "source + diffusion" },
            ],
            legend: "bl",
        }), "fe_masses.png");
    }
}

main();
