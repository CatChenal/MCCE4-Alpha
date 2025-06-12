#!/usr/bin/env python

import logging
from pathlib import Path
from typing import Tuple, Union
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from numpy.polynomial import Polynomial as Poly
import pandas as pd
import warnings
from mcce4.mcce_benchmark import io_utils as iou


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def fit_series(x: pd.Series, y: pd.Series) -> tuple:
    """Fit x to y using Linear Least Square (LLS) fit using
    numpy.Polynomial.fit in 1d.
    Returns:
      A 5-tuple:
       - converged (bool): whether the fit was successful
       - pfit, rmse, r_squared, err_msg
    """
    np.seterr(all="raise")

    converged = True
    pfit = None
    rmse = None
    r_squared = None
    err_msg = None

    try:
        pfit, stats = Poly.fit(x.to_numpy(), y.to_numpy(), 1, full=True)
        # TSS :: Total Sum of Squares
        TSS = np.sum((y - y.mean()) ** 2)
        # rmse: RMS error of the fit
        rmse = np.sqrt(TSS / y.size)
        # r_squared :: coefficent of determination
        # stats[0][0] :: RSS :: Residual Sum of Squares
        r_squared = (1 - stats[0][0] / TSS) + 0.00001

    except (
        IndexError,
        ZeroDivisionError,
        FloatingPointError,
        RuntimeWarning,
        np.AxisError,
    ):
        # e = IndexError('index 0 is out of bounds for axis 0 with size 0')
        # e = FloatingPointError('divide by zero encountered in scalar divide')
        converged = False
        err_msg = "(too few points)"
    except Exception as e:
        # LinAlgError: ("SVD did not converge in Linear Least Squares")
        # RankWarning: Polyfit may be poorly conditioned
        converged = False
        err_msg = str(e)
        if err_msg.startswith("Poly"):
            err_msg = "(bad fit)"
        elif err_msg.startswith("SVD"):
            err_msg = "(unconverged)"
        elif err_msg.startswith("** On entry"):
            err_msg = "(data error, null values?)"

    return converged, pfit, rmse, r_squared, err_msg


def plot_conf_thrup(
    tp_df: pd.DataFrame,
    n_complete: int,
    set_name: str,
    level: int = None,
    out_fp: str = None,
    show: bool = False,
) -> None:
    """Conformers throughput per mcce steps that took over 0 s."""
    # iou.txt2df does not read 1st colname if named index:
    tp_df.rename(columns={"Unnamed: 0": "step", "confs_per_min": "y"}, inplace=True)
    try:
        imax = tp_df.y.idxmax()
    except (AttributeError, ValueError):
        logger.error("No data in confs_per_min column (y)?")
        return

    fig, ax = plt.subplots(1, 1, layout="constrained")
    plt.yscale("log")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylabel("mean conformers/min")
    ax.set_xlabel("mcce steps with dwell time > 0s")

    hdr = f"{set_name.upper()} - Mean conformer throughput per min"
    if level is None:
        hdr = hdr + f"\n(N={n_complete} pdbs)"
    else:
        lev = iou.levels2names[level]
        hdr = hdr + f" @ {lev!r} level\n(N={n_complete} pdbs)"
    fig.suptitle(hdr)

    try:
        markerline, stemlines, baseline = plt.stem(
            tp_df.step, tp_df.y, linefmt="grey", bottom=0
        )
    except AttributeError as e:
        logger.error(f"Plotting failed, missing column?: {e}")
        return

    msize = 4
    plt.setp(
        markerline, ms=msize, markerfacecolor="tab:blue", markeredgecolor="tab:blue"
    )
    plt.setp(stemlines, linewidth=0.5, color="tab:blue")
    plt.setp(baseline, linewidth=0.5, color="k")

    # imax = tp_df.y.idxmax(): set early on
    mx_pt = tp_df.iloc[imax].step, tp_df.iloc[imax].y
    plt.plot(
        tp_df.iloc[imax].step,
        tp_df.iloc[imax].y,
        "o",
        label=f"{mx_pt[1]:,.0f}",
        ms=msize + 1,
        markerfacecolor="tab:red",
        markeredgecolor="tab:red",
    )

    plt.grid(visible=True, which="major", axis="y", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)

    ax.legend(
        title=f"Best: {mx_pt[0].title()}",
        loc="upper center",
        bbox_to_anchor=(0.5, 0.9),
        # borderaxespad=0.0,
        facecolor="tab:cyan",
        framealpha=0.1,
    )

    if out_fp is not None:
        plt.savefig(out_fp)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return


def plot_pkas_fit(
    matched_fp: Path, pks_stats: dict, out_fp: Path = None, show: bool = False,
    comparison: bool=False) -> None:
    """Plot the best fit line of calculated vs experimental or reference pKas
    for the matched pkas in the benchmark.
    Args:
     - matched_fp (Path): path to matched_pkas.csv
     - pks_stats (dict): Output of pkanalysis.matched_pkas_stats(matched_fp)
     - out_fp (Path, None): Filepath to save figure
     - show (bool, False): whether to display the plot
     - comparison (bool, False): For constructing the title
    """
    N = pks_stats["N"]
    matched_df = iou.matched_pks_to_df(matched_fp)
    Y = matched_df.iloc[:, 2]
    X = matched_df.iloc[:, 3]
    # reduced point set for equal & fitted lines:
    eq = np.array([i for i in range(14)])

    unconverged = pks_stats["fit"][0] == False
    if unconverged:
        logger.error(f"The matched pKa data could not be fitted {pks_stats['fit'][1]}.")
        # return
    else:
        if len(pks_stats["fit"]) > 4:
            # a previous version had 5 outputs
            pks_stats["fit"] = pks_stats["fit"][:4]

        b, m, rmse, r_sqr = pks_stats["fit"]
    
    try:
        # a previous version won't have these:
        rmsd, corr = pks_stats["data_stats"]
    except KeyError:
        logger.error("Rerun analysis/comparison: pks_stats is missing rmsd and correlation.",
                     exc_info=1)
        return

    # process df until cols :: res <set1> <set2>
    matched_df["res"] = matched_df.resid.str[:3]
    matched_df.sort_values(by="res", inplace=True)
    matched_df.drop(columns=["pdb", "resid"], inplace=True)
    matched_df = matched_df[["res", matched_df.columns[0], matched_df.columns[1]]]

    # sets names for axes:
    val_cols = matched_df.columns[-2:].tolist()
    res_u = matched_df.res.unique()

    fig, ax = plt.subplots(layout="constrained")
    ax.set_xlim(-1, 14)
    ax.set_ylim(-1, 14)
    which_vals = pks_stats["titr_vals"]

    ax.set_xlabel(f"{X.name.upper()} {which_vals}s; N matched = {N:,}")
    ax.set_ylabel(f"{Y.name.upper()} {which_vals}s")
    ax.spines[["top", "right"]].set_visible(False)

    # the level key from pks_stats is always a tuple; it's a 2-tuple (of tuples)
    # in the context of a comparison.
    levtupl = pks_stats["level"]
    multi_valued = False
    # get title string with level(s) if possible:
    title_str = ""
    if not comparison:
        if levtupl[0]:
            title_str = f"Calculated v. experimental {which_vals} values @ {levtupl[1]!r} level"
        else:
            title_str = f"Calculated v. experimental {which_vals} values"
    else:
        multi_valued = isinstance(levtupl[0], tuple) and len(levtupl[0]) == len(levtupl)
        if multi_valued:
            l1 = levtupl[0][0]
            l2 = levtupl[1][0]
            if not l1 and not l2:
                title_str = f"Compared (Y v. X) {which_vals} values"
            elif l1 == l2:
                title_str = f"Compared (Y v. X) {which_vals} values @ {levtupl[0][1]!r} level"
            else:
                title_str = f"Compared Y (@ {levtupl[0][1]!r}) v. X (@ {levtupl[1][1]!r}) {which_vals} values"
        else:
            title_str = f"Compared (Y v. X) {which_vals} values @ {levtupl[1]!r} level"
            
    fig.suptitle(title_str)

    residues_handles = []
    for name in res_u:
        msk = matched_df.res == name
        x = matched_df[msk][val_cols[1]]
        y = matched_df[msk][val_cols[0]]
        (handle,) = ax.plot(
            x,
            y,
            label=name.upper(),
            marker="o",
            markeredgecolor="black",
            markeredgewidth=.25,
            ls="",
        )
        residues_handles.append(handle)
    
    lines_handles = []
    (equal_line,) = ax.plot(eq, eq, ls="solid", color="r", lw=1.5, label="y = x")

    if not unconverged:
        (fit_line,) = ax.plot(
            eq,
            m * eq + b,
            label=r"$"
            + f"y = {m:.2f}x {b:+.2f}, R^2: {r_sqr:.2f}$ \n $RMSD: {rmsd:.2f}, Corr: {corr:.2f}$",
            ls="-",
            lw=1.5,
            color="k",
        )
    else:
        # dummy fit line
        (fit_line,) = ax.plot(
            0, 0,
            label= f"No fit {pks_stats['fit'][1]}\nRMSD: {rmsd:.2f}, Corr: {corr:.2f}",
            ls="-",
            lw=1.5,
            color="k",
        )
    
    lines_handles.extend([fit_line, equal_line])

    cm = plt.get_cmap("tab20")
    # bounds = [3*titr_d, 2*titr_d, 1.5*titr_d, titr_d, 0.5*titr_d]
    for c, v in enumerate(pks_stats["bounds"]):
        if c in [2, 4]:
            # fractional bounds, don't show
            continue
        (bound_handle1,) = ax.plot(eq, eq + v, ls=":", color=cm(c))
        ax.plot(eq, eq - v, ls=":", color=cm(c))
        # Calc %points within these bounds around equality:
        points_within = np.logical_and(Y >= X - v, Y <= X + v)
        percentage = np.sum(points_within) / N
        # Set the label for each upper bound handle
        lblstr = r"$\pm$" + f"{v:.0f}: {percentage:> 8.2%}"
        bound_handle1.set_label(lblstr)
        lines_handles.append(bound_handle1)

    # Add a separate legend for lines:
    ax.add_artist(
        ax.legend(
            handles=lines_handles,
            frameon=False,
            loc="upper left",
        )
    )
    lc = 6
    if res_u.size > 6:
        lc = 5
    ax.legend(
        handles=residues_handles, ncol=lc,
        frameon=False,
        bbox_to_anchor=(0.5, -0.13), loc="upper center",
    )

    if out_fp is not None:
        plt.savefig(out_fp)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return


def correl(x: pd.Series, y: pd.Series, somename: str = None) -> float:
    """Compute the correlation between x and y, two pandas.Series.
    Logs a warning if correlation reset to 0.
    """
    if somename is None:
        msg1 = "Correlation set to 0."
    else:
        msg1 = f"Correlation for {somename} set to 0."
    
    if y.size <= 1:
        logger.warning(msg1 + " Too few points.")
        return 0.
    
    if all(abs(y - x) == 0):
        # Same values; no change
        return 1.

    np.seterr(all="raise")
    with warnings.catch_warnings():
        warnings.filterwarnings("error")
        try:
            corr = y.corr(x)
        except (RuntimeWarning, FloatingPointError) as e:
            corr = 0
            logger.warning(msg1)

    return corr


def plot_res_analysis(
    matched_fp: Path,
    res_stats: dict,
    level: Union[int, Tuple[int]] = None,
    out_fp: str = None,
    show: bool = False,
    figsize: tuple = (9, 11),
) -> None:
    """
    Plot the best fit line of matched pKas csv file grouped by residue.
    matched_fp (Path): path to matched_pkas.csv; columns: pdb resid <set1> <set2> delta.
    In the context of a comparison, the level argument will be a 2-tuple if the conformer-
    making level in the two benchmark folders differ.
    """
    matched_df = iou.matched_pks_to_df(matched_fp)
    # A, B = matched_df.columns[2:4]  # A :: calc; B :: ref or expl

    # process df until cols :: res <set1> <set2>
    matched_df["res"] = matched_df.resid.str[:3]
    matched_df.sort_values(by="res", inplace=True)
    matched_df.drop(columns=["pdb", "resid"], inplace=True)
    matched_df = matched_df[["res", matched_df.columns[0], matched_df.columns[1]]]

    # sets names for axes:
    val_cols = matched_df.columns[-2:].tolist()
    res_u = matched_df.res.unique()
    N_res = len(res_u)
    num_cols = 2
    num_rows = int(np.ceil(N_res / num_cols))

    fig = plt.figure(figsize=figsize, layout="constrained")
    supti = "Residue Analysis"
    if level is None:
        fig.suptitle(supti)
    else:
        if isinstance(level, tuple):
            lname1 = iou.levels2names.get(level[0])
            if lname1 is None:
                lname1 = "unknown"
            lname2 = iou.levels2names.get(level[1])
            if lname2 is None:
                lname2 = "unknown"
            fig.suptitle(supti + f" @ {lname1!r}|{lname2!r} levels")
        else:
            lname1 = iou.levels2names.get(level) or "unknown"
            #if lname1 is None:
            #    lname1 = "unknown"
            fig.suptitle(supti + f" @ {lname1!r} level")

    cm = plt.get_cmap("tab20")
    gs = plt.GridSpec(num_rows, num_cols, figure=fig)

    for i, name in enumerate(res_u):
        r = i // num_cols
        c = i % num_cols
        ax = fig.add_subplot(gs[r, c])
        msk = matched_df.res == name
        x = matched_df[msk][val_cols[1]]
        y = matched_df[msk][val_cols[0]]
        
        corr = correl(x, y, name)
        rmsd = np.sqrt(np.mean((y - x)**2))

        ax.plot(x, x, ls="solid", lw=1.5, color="k", label="y = x")

        converged, pfit, rmse, r_sqr, err_msg = fit_series(x, y)
        if converged:
            b, m = pfit.convert()
            ax.plot(
                x,
                pfit(x),
                "-",
                color="r",
                lw=1.5,
                ms=2,
                label=r"$"
                + f"y = {m:.2f}x {b:+.2f}, R^2: {r_sqr:.2f}$ \n $RMSD: {rmsd:.2f}, Corr: {corr:.2f}$",
            )
            ax.set_title(f"{name.upper()}", y=0.95, pad=-12)
        else:
            ax.set_title(f"{name}: {err_msg}\nRMSD: {rmsd:.2f}, Corr: {corr:.2f}",
                          y=0.95, pad=-12)

        # get the res stdev from the dict:
        yerrs = [res_stats[name]["delta_std"]] * len(y)
        plt.errorbar(
            x, y, yerr=yerrs,
            marker="o", 
            markeredgecolor="k",
            markeredgewidth=.3,
            ls="", color=cm(i + 1),
            #label=name
        )

        xlbl = val_cols[1] if (r == num_rows - 1) else ""
        ax.xaxis.set(label_text=xlbl)
        ylbl = val_cols[0] if (c == 0) else ""
        ax.yaxis.set(label_text=ylbl)

        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(framealpha=0.5,
                  loc="lower center",
                  fontsize="small",
                  )

    if out_fp is not None:
        plt.savefig(out_fp)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return


def plot_sumcrg_diff(
    sc_diff12: pd.DataFrame, pdbid: str, sets_names: tuple, out_fp: str = None, show: bool = False
):
    """Plot the sum_crg difference at a titration point for a 'pdbid'.
    Args:
      - sc_diff12 (pd.DataFrame): the output of `comparison.get_allsumcrg_diff(sc1, sc2,
        titr_point=x)` where sc1 and sc2 are the collated sum_crg files of the two compared
        sets of runs.
      - pdbid (str): the pdbid used to filter sc_diff12.
      - sets_names (tuple): the name of the data sources.
      - out_fp (str): the file path to save the figure.
      - show (bool): whether to display the plot.
    WIP
    prot_df = all_sumcrg_df[all_sumcrg_df.Diff == pdbid]
    msk = abs(prot_df[prot_df.columns[-1]]) >= scdiff_thresh
    if not prot_df[msk].size:
        print(f"No significant sum_crg difference in {pdbid}.")
        continue

    """
    df = sc_diff12[sc_diff12.Diff == pdbid]
    if not df.size:
        logger.error(f"This pdbid: {pdbid} was not found.")
        return

    conf_col = df.columns[1]
    delta_col = df.columns[2]
    df.loc[:, conf_col] = df[conf_col].str.replace(r"_$", "", regex=True)

    totals = ["Net_Charge", "Protons", "Electrons"]
    tot_msk = df[conf_col].isin(totals)
    totals_vals = df.loc[tot_msk, delta_col].values

    df = df.drop(df[tot_msk].index, axis=0)
    N = df.shape[0]

    offs = 0.0
    if N < 63:
        offs = 4.0
    H = np.ceil(N * 10 / 63 + offs)
    # figsize=(5, H)
    fig, ax = plt.subplots(layout="constrained")
    plt.grid(alpha=0.2)

    colors = ["tab:red", "grey", "tab:blue"]
    colr = np.where(
        df[delta_col] < 0, colors[0], np.where(df[delta_col] > 0, colors[2], colors[1])
    )
    plt.scatter(
        df[delta_col],
        df[conf_col],
        marker="o",
        s=10,
        c=colr,
        label=pdbid,
    )

    plt.axvline(color="grey", ymax=0.98, linestyle="--", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel("Charge difference")
    ax.xaxis.set_major_formatter("{x:.2f}")
    ax.set_title(f"Net charge: {sets_names[0]} - {sets_names[1]} @ {conf_col} {delta_col} | {pdbid}",
                 y=1.0,
                 # fontsize="medium",
                 )

    leg_handles = [
        Line2D([], [], ls="", marker="", color="k", label=f"{totals[0]}: {totals_vals[0]:.1f}"),
        Line2D([], [], ls="", marker="", color="tab:blue", label=f"{totals[1]}: {totals_vals[1]:.1f}"),
        Line2D([], [], ls="", marker="", color="tab:pink", label=f"{totals[2]}: {totals_vals[2]:.1f}")
    ]
    ax.legend(
        handles=leg_handles,
        handlelength=.5,
        handletextpad=.5,
        bbox_to_anchor=(0.5, -0.1),
        loc="upper center",
        ncol=3,
        fontsize="medium",
        markerfirst=False,
        framealpha=1,
        edgecolor="white",
        borderaxespad=0.5,
        labelcolor="linecolor",
    )

    if out_fp is not None:
        plt.savefig(out_fp)
    if show:
        plt.show()
    else:
        plt.close(fig)

    return
