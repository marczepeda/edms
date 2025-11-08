import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm, truncnorm
import argparse
from rich_argparse import RichHelpFormatter
from rich import print as rprint

def truncnorm_pdf(x, mu, sigma, lo=0.0, hi=1.0):
    a = (lo - mu) / sigma
    b = (hi - mu) / sigma
    return truncnorm(a, b, loc=mu, scale=sigma).pdf(x)

def main():
    # Custom formatter for rich help messages
    class MyFormatter(RichHelpFormatter):
        styles = {
            "argparse.prog": "green",           # program name
            "argparse.args": "cyan",            # positional arguments
            "argparse.option": "",              # options like --flag
            "argparse.metavar": "dark_magenta", # meta variable (actual function argument name)
            "argparse.help": "blue",            # help text
            "argparse.text": "green",           # normal text in help message
            "argparse.groups": "red",           # group titles
            "argparse.description": "",         # description at the top
            "argparse.epilog": "",              # ... -h; epilog at the bottom
            "argparse.syntax": "white",         # []
        }

    parser = argparse.ArgumentParser(description="Calculate and plot mutation frequency distributions.",
                                     epilog="Example usage: python rate.py --n-mut 4000 --mean 0.20 --sd 0.02 --qphred 25",
                                     formatter_class=MyFormatter)
    parser.add_argument('--n-mut', type=int, nargs='+', default=[4000], help='Number(s) of designed mutations; can pass multiple values (default: 4000)')
    parser.add_argument('--mean', type=float, nargs='+', default=[0.20], help='Mean total edit rate(s) (fraction), space-separated (default: 0.20)')
    parser.add_argument('--sd', type=float, default=0.02, help='Standard deviation of total edit rate (fraction) (default: 0.02)')
    parser.add_argument('--qphred', type=int, nargs='+', default=[25], help='Phred quality score(s); can specify multiple values, e.g. --qphred 25 30 35 (default: 25)')
    parser.add_argument('--min-count', type=int, default=10, help='Target counts per average individual mutation (x). Computes reads required so expected counts ≥ x (default: 10)')
    parser.add_argument('--amplicon-len', type=int, nargs='+', default=None,
                        help='Amplicon length(s) in bases; if provided, estimate reads lost to expected per-read errors at each Q score')
    parser.add_argument('--no-plot', action='store_true', help='Skip plotting for headless runs')
    parser.add_argument('--show', action='store_true', help='Show plot interactively (default: False)', default=False)
    parser.add_argument('--save', action='store_true', help='Save output files', default=False)
    args = parser.parse_args()

    N_MUTS = args.n_mut          # list of numbers of designed mutations
    MEANS = args.mean            # list of mean total edit rates (fractions)
    SD   = args.sd           # sd of total edit rate (fraction)
    Q_PHREDS = args.qphred
    TARGET_COUNT = args.min_count
    AMPLICON_LENS = args.amplicon_len

    # Per-base error probability for Phred+33 = Q_PHRED
    p_errs = [10 ** (-q / 10.0) for q in Q_PHREDS]

    for MU in MEANS:
        for N_MUT in N_MUTS:
            mean_F = MU / N_MUT
            sd_F   = SD / N_MUT
            q_lo_R = float(np.clip(norm.ppf(0.025, loc=MU, scale=SD), 0.0, 1.0))
            q_hi_R = float(np.clip(norm.ppf(0.975, loc=MU, scale=SD), 0.0, 1.0))
            q_lo_F = q_lo_R / N_MUT
            q_hi_F = q_hi_R / N_MUT
            print(f"[μ={MU:.2%}, N={N_MUT}] Single mutation mean ≈ {mean_F*100:.6f}%  sd ≈ {sd_F*100:.6f}%")
            print(f"[μ={MU:.2%}, N={N_MUT}] Single mutation 95% ≈ [{q_lo_F*100:.6f}%, {q_hi_F*100:.6f}%]")

    for q, p in zip(Q_PHREDS, p_errs):
        print(f"Per-base error at Q{q} ≈ {p*100:.6f}%")

    # Compute reads required to observe at least TARGET_COUNT counts for the average individual mutation
    reads_rows = []
    print(f"\nReads required to observe ≥ {TARGET_COUNT} counts for the average single mutation:")
    for MU in MEANS:
        for N_MUT in N_MUTS:
            p_single = MU / N_MUT  # average single-mutation frequency (fraction)
            if p_single <= 0:
                reads_required = np.inf
            else:
                reads_required = int(np.ceil(TARGET_COUNT / p_single))
            reads_rows.append({
                'mu_total': MU,
                'N_mut': N_MUT,
                'target_count': TARGET_COUNT,
                'p_single': p_single,
                'reads_required': reads_required
            })
            print(f"[μ={MU:.2%}, N={N_MUT}] p_single≈{p_single:.6e} ⇒ reads≈{reads_required:,}")

    # If amplicon lengths are provided, estimate reads lost due to per-read error
    loss_rows = []
    if AMPLICON_LENS:
        print("\nExpected reads lost to sequencing errors (per-read) given amplicon length(s):")
        for MU in MEANS:
            for N_MUT in N_MUTS:
                # recompute reads_required for this (MU, N_MUT)
                p_single = MU / N_MUT
                reads_required = np.inf if p_single <= 0 else int(np.ceil(TARGET_COUNT / p_single))
                for L in AMPLICON_LENS:
                    for q, p_base in zip(Q_PHREDS, p_errs):
                        # Probability a read of length L has ≥1 error (assuming independence)
                        p_read_error = 1.0 - (1.0 - p_base) ** L
                        lost = int(np.round(reads_required * p_read_error)) if np.isfinite(reads_required) else np.inf
                        effective = reads_required - lost if np.isfinite(reads_required) else np.inf
                        loss_rows.append({
                            'mu_total': MU,
                            'N_mut': N_MUT,
                            'target_count': TARGET_COUNT,
                            'amplicon_len': L,
                            'qphred': q,
                            'p_base_error': p_base,
                            'p_read_error': p_read_error,
                            'reads_required': reads_required,
                            'expected_lost_reads': lost,
                            'expected_effective_reads': effective
                        })
                        print(f"[μ={MU:.2%}, N={N_MUT}, L={L}bp, Q{q}] p_read_error≈{p_read_error:.6%} ⇒ lost≈{lost:,} (effective≈{effective:,})")

    if not args.no_plot:
        nrows = len(N_MUTS)
        ncols = len(MEANS)
        fig_w = 6 * ncols
        fig_h = 4.5 * nrows
        fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h), squeeze=False)

        for i, N_MUT in enumerate(N_MUTS):
            for j, MU in enumerate(MEANS):
                ax = axes[i, j]
                # Compute the per-mutation PDF for this (N, MU)
                x_F = np.linspace(0, 1.0 / N_MUT, 2000)
                pdf_F = truncnorm_pdf(N_MUT * x_F, MU, SD, 0.0, 1.0) * N_MUT
                ax.plot(x_F * 100, pdf_F, linewidth=2)
                # Add color fill under the PDF curve
                ax.fill_between(x_F * 100, pdf_F, color='skyblue', alpha=0.3)

                # Vertical reference lines: per-base error at each Q with inline labels
                for q, p in zip(Q_PHREDS, p_errs):
                    x_pos = p * 100
                    ax.axvline(x_pos, linestyle="--", linewidth=1.2)
                    # Place a rotated Q label near the top of the axis, aligned to the line
                    y_top = ax.get_ylim()[1]
                    ax.text(x_pos, y_top * 0.98, f"Q{q}", rotation=90,
                            va="top", ha="center", fontsize=8,
                            bbox=dict(facecolor="white", alpha=0.6, edgecolor="none"))

                # Axes labels and titles
                if i == nrows - 1:
                    ax.set_xlabel("Rate (%)")
                if i == 0:
                    ax.set_title(f"Aggregate Mutation Rate = {MU:.1%}")
                if j == 0:
                    ax.set_ylabel(f"N = {N_MUT}\nDensity (F)")
                ax.grid(True, which="both", linestyle=":", linewidth=0.5)


        fig.suptitle("Single-mutation frequency vs. per-base error rate", fontsize=14, y=0.99)
        plt.tight_layout()

        # Second figure: Reads required vs reads lost (requires amplicon lengths)
        loss_plot_created = False
        fig2 = None
        if AMPLICON_LENS:
            df_loss = pd.DataFrame(loss_rows)
            if not df_loss.empty:
                nrows2 = len(N_MUTS)
                ncols2 = len(AMPLICON_LENS)
                fig2_w = 7 * ncols2
                fig2_h = 3 * nrows2
                fig2, axes2 = plt.subplots(nrows2, ncols2, figsize=(fig2_w, fig2_h), squeeze=False)
                for i, N_MUT in enumerate(N_MUTS):
                    for j, L in enumerate(AMPLICON_LENS):
                        ax2 = axes2[i, j]
                        sub = df_loss[(df_loss['N_mut'] == N_MUT) & (df_loss['amplicon_len'] == L)]
                        if sub.empty:
                            ax2.set_visible(False)
                            continue
                        # Stacked horizontal bars per Q, grouped by mean (μ)
                        qs = sorted(sub['qphred'].unique())
                        mus = sorted(sub['mu_total'].unique())
                        if not qs or not mus:
                            ax2.set_visible(False)
                            continue
                        y_base = np.arange(len(qs))
                        bar_height = 0.8 / max(1, len(mus))
                        offsets = np.linspace(-0.4 + bar_height/2, 0.4 - bar_height/2, len(mus))
                        xmax = 0
                        for k, MU_ in enumerate(mus):
                            subM = sub[sub['mu_total'] == MU_]
                            # Ensure order by Q
                            req = []
                            lost = []
                            for q in qs:
                                row = subM[subM['qphred'] == q]
                                if row.empty:
                                    req.append(0)
                                    lost.append(0)
                                else:
                                    rreq = float(row['reads_required'].iloc[0])
                                    rlost = float(row['expected_lost_reads'].iloc[0])
                                    req.append(rreq)
                                    lost.append(rlost)
                            y_pos = y_base + offsets[k]
                            ax2.barh(y_pos, req, height=bar_height)
                            ax2.barh(y_pos, lost, height=bar_height, left=req)
                            # Annotate total at the end of each stacked bar and tag μ
                            totals = [a + b for a, b in zip(req, lost)]
                            for yv, tot, a, b in zip(y_pos, totals, req, lost):
                                if tot > 0:
                                    ax2.text(tot * 1.01, yv, f"{tot/10**6:.1f}M", va='center', fontsize=10)
                                if a > 0:
                                    ax2.text(a * 0.5, yv, f"required (μ={MU_:.0%})", va='center', ha='center', fontsize=10)
                                if b > 0:
                                    ax2.text(a + b * 0.5, yv, "lost", va='center', ha='center', fontsize=10)
                            xmax = max(xmax, max(totals) if totals else 0)
                        # Y ticks at Q categories
                        ax2.set_yticks(y_base)
                        ax2.set_yticklabels([f"Q{q}" for q in qs])
                        # Labels
                        if i == nrows2 - 1:
                            ax2.set_xlabel('Reads (required + lost)')
                        if j == 0:
                            ax2.set_ylabel('Phred+33 score')
                        # Title now includes amplicon length
                        ax2.set_title(f"L={L} bp | N={N_MUT} | n={TARGET_COUNT}")
                        ax2.set_xlim(0, xmax * 1.1 if xmax > 0 else 1)
                        ax2.grid(True, axis='x', linestyle=':', linewidth=0.5)
                fig2.suptitle('Reads required + lost to errors for PE screens', fontsize=14, y=0.98)
                plt.tight_layout()
                loss_plot_created = True

    if args.show:
        plt.show()

    if args.save:
        out_png = "mutation_frequency_distribution_facets.png"
        fig.savefig(out_png, dpi=300)
        rprint(f"[green]Wrote:[/green] {out_png}")
        # Save second figure if present
        try:
            if 'fig2' in locals() and fig2 is not None:
                out_png2 = "reads_required_vs_reads_lost.png"
                fig2.savefig(out_png2, dpi=300)
                rprint(f"[green]Wrote:[/green] {out_png2}")
        except Exception:
            pass

    # Always save per-combination CSVs when --save is set
    if args.save:
        wrote = []
        # Save a total-rate PDF for each mean
        x_R = np.linspace(0, 1, 2000)
        for MU in MEANS:
            pdf_R = truncnorm_pdf(x_R, MU, SD, 0.0, 1.0)
            out_total = f"total_rate_pdf_mu{MU:.2f}.csv"
            pd.DataFrame({"rate_total_percent": x_R * 100, "pdf_total": pdf_R}).to_csv(out_total, index=False)
            wrote.append(out_total)
        # Save single-mutation PDFs for each (mean, N) pair
        for MU in MEANS:
            for N_MUT in N_MUTS:
                x_F = np.linspace(0, 1.0 / N_MUT, 2000)
                pdf_F = truncnorm_pdf(N_MUT * x_F, MU, SD, 0.0, 1.0) * N_MUT
                out_single = f"single_mutation_pdf_mu{MU:.2f}_N{N_MUT}.csv"
                pd.DataFrame({"rate_single_mut_percent": x_F * 100, "pdf_single_mut": pdf_F}).to_csv(out_single, index=False)
                wrote.append(out_single)
        # Save reads required table for target counts per single mutation
        out_reads = f"reads_required_target{TARGET_COUNT}.csv"
        pd.DataFrame(reads_rows).to_csv(out_reads, index=False)
        wrote.append(out_reads)
        if AMPLICON_LENS:
            out_loss = f"reads_lost_due_to_errors_target{TARGET_COUNT}.csv"
            pd.DataFrame(loss_rows).to_csv(out_loss, index=False)
            wrote.append(out_loss)
        rprint("[green]Wrote:[/green] " + ", ".join(wrote))

if __name__ == "__main__":
    main()