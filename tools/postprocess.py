#!/usr/bin/env python3
"""
tools/postprocess.py - Post-processing of control_print_tabular output.

Reads the SQLite database produced by tochnog's control_print_tabular
(or a plain CSV) and provides:

  info   - summary of the database (tables, dofs, time steps, nodes)
  stats  - descriptive statistics of one or more dofs
  plot   - time series plot of a dof (one or more nodes, saved to PNG)
  line   - profile of a dof along a geometry_line (interpolated)
  user   - compute a Python variable per (node,t) and store it in user_data

Requires: pandas, numpy, matplotlib (plot/line only need matplotlib).
Usage examples:
  postprocess.py info biax20.sqlite
  postprocess.py stats biax20.sqlite --dof sigyy_11
  postprocess.py plot biax20.sqlite --dof sigyy_11 --node 0 --out sig.png
  postprocess.py line biax20.sqlite --dof sigyy_11 --x1 0 --y1 0 --x2 1 --y2 0
  postprocess.py user biax20.sqlite --name kin --expr "0.5*(velx_0**2+vely_0**2)"
"""

import argparse
import sqlite3
import sys

try:
    import pandas as pd
except ImportError:
    pd = None


# ---------------------------------------------------------------- loading


def load_sqlite(path):
    """Return (frames, meta) where frames is a dict name->DataFrame."""
    con = sqlite3.connect(path)
    meta = {}
    try:
        for k, v in con.execute("SELECT key,value FROM meta"):
            meta[k] = v
    except sqlite3.Error:
        pass

    frames = {}
    try:
        frames["primary"] = pd.read_sql("SELECT * FROM primary_data", con)
    except sqlite3.Error:
        pass
    for table in ("derived", "coords", "user_data"):
        try:
            frames[table] = pd.read_sql(f"SELECT * FROM {table}", con)
        except sqlite3.Error:
            frames[table] = None
    con.close()
    return frames, meta


def pivot_primary(prim):
    """Long primary_data -> wide (node, t) x dof columns."""
    return prim.pivot_table(
        index=["node", "t"], columns="dof", values="value", aggfunc="first"
    ).reset_index()


def wide_with_derived(frames):
    """Pivot primary and merge derived columns (same (node,t) key)."""
    wide = pivot_primary(frames["primary"])
    der = frames.get("derived")
    if der is not None and len(der):
        wide = wide.merge(der, on=["node", "t"], how="left")
    return wide


def require_pandas():
    if pd is None:
        sys.exit("Error: pandas is required. Install with: pip install pandas")


# ---------------------------------------------------------------- subcommands


def cmd_info(args):
    frames, meta = load_sqlite(args.db)
    print(f"Database: {args.db}")
    print(f"meta: {meta}")
    for name, df in frames.items():
        if df is None:
            continue
        print(f"\ntable {name}: {len(df)} rows")
        if name == "primary":
            dofs = sorted(df["dof"].unique())
            print(f"  dofs ({len(dofs)}): {', '.join(map(str, dofs))}")
            print(
                f"  times: {len(df['t'].unique())} steps, "
                f"t in [{df['t'].min():g}, {df['t'].max():g}]"
            )
            print(f"  nodes: {sorted(df['node'].unique())}")
        elif name == "derived":
            print(f"  columns: {', '.join(df.columns)}")
        elif name == "coords":
            print(f"  nodes: {len(df)}")
        elif name == "user_data" and len(df):
            print(f"  variables: {sorted(df['variable'].unique())}")


def cmd_stats(args):
    require_pandas()
    frames, _ = load_sqlite(args.db)
    wide = wide_with_derived(frames)

    if args.dof:
        cols = [c for c in args.dof if c in wide.columns]
        missing = [c for c in args.dof if c not in wide.columns]
        if missing:
            print(f"Warning: dofs not found: {missing}")
        if not cols:
            sys.exit("Error: none of the requested dofs found.")
    else:
        cols = [c for c in wide.columns if c not in ("node", "t")]
        if not cols:
            sys.exit("Error: no dof columns to compute stats on.")

    stats = wide[cols].describe().T
    print(stats.to_string())
    if args.csv:
        stats.to_csv(args.csv)
        print(f"\nstats written to {args.csv}")


def cmd_plot(args):
    require_pandas()
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    frames, meta = load_sqlite(args.db)
    wide = wide_with_derived(frames)

    if args.dof not in wide.columns:
        sys.exit(
            f"Error: dof '{args.dof}' not found. Available: "
            f"{[c for c in wide.columns if c not in ('node', 't')]}"
        )

    nodes = args.node if args.node else sorted(wide["node"].unique())
    if args.node:
        nodes = [n for n in nodes if n in set(wide["node"])]
        if not nodes:
            sys.exit("Error: none of the requested nodes found.")

    plt.figure(figsize=(8, 5))
    for n in nodes:
        sub = wide[wide["node"] == n].sort_values("t")
        plt.plot(sub["t"], sub[args.dof], marker=".", label=f"node {n}")
    plt.xlabel("t")
    plt.ylabel(args.dof)
    plt.title(f"{args.dof} vs t ({args.db})")
    plt.legend()
    plt.grid(True, alpha=0.3)
    out = args.out or f"{args.dof}.png"
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    print(f"plot saved to {out}")


def cmd_line(args):
    require_pandas()
    import numpy as np

    frames, _ = load_sqlite(args.db)
    coords = frames.get("coords")
    if frames.get("primary") is None or not len(frames["primary"]):
        sys.exit("Error: no primary_data in database.")
    if coords is None or not len(coords):
        sys.exit("Error: no coords table in database.")
    wide = wide_with_derived(frames)
    wide = wide.merge(coords, on="node", how="left")

    if args.dof not in wide.columns:
        sys.exit(f"Error: dof '{args.dof}' not found.")

    t = args.t
    if t is None:
        t = wide["t"].max()
    if t not in set(wide["t"]):
        sys.exit(
            f"Error: time {t} not present. Available times near: "
            f"{sorted(set(wide['t']))[:8]}"
        )

    sub = wide[wide["t"] == t].copy()
    sub = sub[sub[args.dof].notna()]

    # geometry line from p0 to p1 (only x,y used)
    p0 = np.array([args.x1, args.y1])
    p1 = np.array([args.x2, args.y2])
    d = p1 - p0
    length = np.linalg.norm(d)
    if length == 0:
        sys.exit("Error: zero-length geometry line.")
    u = d / length

    xy = sub[["x", "y"]].to_numpy()
    # projection of each node onto the line (0..1 parameter)
    s = ((xy - p0) @ u) / length
    # perpendicular distance
    perp = np.linalg.norm(xy - (p0 + np.outer(s, u) * length), axis=1)

    tol = args.tol
    sel = (s >= 0) & (s <= 1) & (perp <= tol)
    if sel.sum() == 0:
        sys.exit(
            "Error: no nodes within the tolerance of the line. "
            "Increase --tol or check the line coordinates."
        )
    res = sub.loc[sel, ["node", "x", "y", args.dof]].copy()
    res["s"] = s[sel] * length  # distance along the line
    res = res.sort_values("s")

    if args.csv:
        res.to_csv(args.csv, index=False)
        print(f"profile written to {args.csv}")

    print(res.to_string(index=False))
    if args.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 5))
        plt.plot(res["s"], res[args.dof], marker="o")
        plt.xlabel("distance along line")
        plt.ylabel(args.dof)
        plt.title(f"{args.dof} along line at t={t}")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(args.plot, dpi=150)
        print(f"plot saved to {args.plot}")


def cmd_user(args):
    require_pandas()
    frames, _ = load_sqlite(args.db)
    if frames.get("primary") is None or not len(frames["primary"]):
        sys.exit("Error: no primary_data in database.")
    wide = wide_with_derived(frames)

    expr = args.expr
    cols = [c for c in wide.columns if c not in ("node", "t")]
    # validate the expression references existing columns only
    import re

    tokens = set(re.findall(r"[A-Za-z_]\w*", expr))
    unknown = (
        tokens
        - set(cols)
        - {
            "node",
            "t",
            "abs",
            "sqrt",
            "exp",
            "log",
            "sin",
            "cos",
            "tan",
            "min",
            "max",
            "sum",
            "mean",
            "len",
        }
    )
    if unknown:
        sys.exit(f"Error: expression references unknown columns: {unknown}")
    try:
        wide[args.name] = wide.eval(expr)
    except Exception as e:
        sys.exit(f"Error evaluating expression: {e}")

    long = wide[["node", "t", args.name]].melt(
        id_vars=["node", "t"],
        value_vars=[args.name],
        var_name="variable",
        value_name="value",
    )
    long = long[long["value"].notna()]
    long = long[["node", "variable", "t", "value"]]

    con = sqlite3.connect(args.db)
    cur = con.cursor()
    cur.execute(
        "CREATE TABLE IF NOT EXISTS user_data ("
        " node INTEGER, variable TEXT, t REAL, value REAL,"
        " PRIMARY KEY (node, variable, t));"
    )
    rows = [tuple(r) for r in long.itertuples(index=False)]
    cur.executemany(
        "INSERT OR REPLACE INTO user_data (node,variable,t,value) VALUES (?,?,?,?)",
        rows,
    )
    con.commit()
    n = len(rows)
    con.close()
    print(f"stored {n} user_data rows (variable '{args.name}').")


# ---------------------------------------------------------------- main


def main():
    p = argparse.ArgumentParser(
        description="Post-processing of tochnog control_print_tabular output."
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    def add_db(sp):
        sp.add_argument("db", help="path to the .sqlite database")

    sp = sub.add_parser("info", help="summary of the database")
    add_db(sp)

    sp = sub.add_parser("stats", help="descriptive statistics of dofs")
    add_db(sp)
    sp.add_argument(
        "--dof", action="append", help="dof to describe (repeatable); default: all"
    )
    sp.add_argument("--csv", help="also write stats to a CSV file")

    sp = sub.add_parser("plot", help="time series plot of a dof")
    add_db(sp)
    sp.add_argument("--dof", required=True, help="dof to plot")
    sp.add_argument(
        "--node", type=int, action="append", help="node(s) to plot; default: all nodes"
    )
    sp.add_argument("--out", help="output PNG path")

    sp = sub.add_parser("line", help="dof profile along a geometry line")
    add_db(sp)
    sp.add_argument("--dof", required=True, help="dof to profile")
    sp.add_argument("--x1", type=float, required=True)
    sp.add_argument("--y1", type=float, required=True)
    sp.add_argument("--x2", type=float, required=True)
    sp.add_argument("--y2", type=float, required=True)
    sp.add_argument("--t", type=float, help="time step (default: last)")
    sp.add_argument(
        "--tol",
        type=float,
        default=0.1,
        help="perpendicular tolerance to select nodes (default 0.1)",
    )
    sp.add_argument("--csv", help="write profile to a CSV file")
    sp.add_argument("--plot", help="save profile plot to PNG")

    sp = sub.add_parser("user", help="compute a variable and store in user_data")
    add_db(sp)
    sp.add_argument("--name", required=True, help="variable name (column in user_data)")
    sp.add_argument(
        "--expr",
        required=True,
        help="Python expression over dof columns, e.g. "
        "'sigyy_11 + 200' or '0.5*(velx_0**2+vely_0**2)'",
    )

    args = p.parse_args()
    fn = {
        "info": cmd_info,
        "stats": cmd_stats,
        "plot": cmd_plot,
        "line": cmd_line,
        "user": cmd_user,
    }[args.cmd]
    fn(args)


if __name__ == "__main__":
    main()
