#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
import sys
import argparse

import openmm as mm
from openmm import unit, Platform
from openmm.app import (
    AmberPrmtopFile, AmberInpcrdFile,
    Simulation, DCDReporter, StateDataReporter, PDBFile,
    NoCutoff, HBonds, OBC2
)

def pick_platform(requested: str | None):
    """
    Try to pick a platform. If requested is given, use it.
    Otherwise prefer GPU platforms if available, fall back to CPU.
    """
    if requested:
        p = Platform.getPlatformByName(requested)
        props = {}
        if requested in ("CUDA", "OpenCL", "HIP", "Metal"):
            props["Precision"] = "mixed"
        return p, props

    for name in ("CUDA", "HIP", "OpenCL", "Metal", "CPU"):
        try:
            p = Platform.getPlatformByName(name)
            props = {}
            if name in ("CUDA", "HIP", "OpenCL", "Metal"):
                props["Precision"] = "mixed"
            return p, props
        except Exception:
            continue

    # final fallback
    return Platform.getPlatformByName("CPU"), {}

def ns_to_steps(ns: float, dt_fs: float) -> int:
    dt_ps = dt_fs / 1000.0
    total_ps = ns * 1000.0
    return int(round(total_ps / dt_ps))

def ps_to_steps(ps: float, dt_fs: float) -> int:
    dt_ps = dt_fs / 1000.0
    return int(round(ps / dt_ps))

def main():
    ap = argparse.ArgumentParser(description="OpenMM GB (implicit solvent) MD from Amber prmtop/inpcrd")
    ap.add_argument("--prmtop", default="system.prmtop", help="Amber topology file")
    ap.add_argument("--inpcrd", default="system.inpcrd", help="Amber coordinate/restart file")

    ap.add_argument("--out_pdb", default="ref.pdb", help="Output reference PDB (after heating, before production)")
    ap.add_argument("--out_dcd", default="traj.dcd", help="Output DCD trajectory")
    ap.add_argument("--time_txt", default="time.txt", help="Total wall time file")

    ap.add_argument("--dt_fs", type=float, default=2.0, help="Timestep (fs), default 2.0")
    ap.add_argument("--friction", type=float, default=1.0, help="Langevin friction (1/ps), default 1.0")
    ap.add_argument("--seed", type=int, default=2025, help="Random seed")

    ap.add_argument("--min_iter", type=int, default=2000, help="Minimization max iterations")
    ap.add_argument("--min_tol_kj_per_mol_nm", type=float, default=10.0,
                    help="Minimization tolerance (kJ/mol/nm), default 10")

    ap.add_argument("--t_start", type=float, default=10.0, help="Heating start temperature (K), default 10")
    ap.add_argument("--t_target", type=float, default=340.0, help="Target temperature (K), default 330")
    ap.add_argument("--heat_ps", type=float, default=340.0, help="Heating duration (ps), default 330 ps")

    ap.add_argument("--prod_ns", type=float, default=0.5, help="Production duration (ns), default 0.5 ns")
    ap.add_argument("--frame_interval_ns", type=float, default=0.005,
                    help="Trajectory frame interval (ns), default 0.005 ns (=5 ps)")

    ap.add_argument("--salt_molar", type=float, default=0.15,
                    help="Implicit solvent salt concentration (M), default 0.15 (physiological)")

    ap.add_argument("--platform", default=None,
                    help="Force a platform name: CUDA/OpenCL/HIP/Metal/CPU (optional)")

    args = ap.parse_args()

    prmtop = AmberPrmtopFile(args.prmtop)
    inpcrd = AmberInpcrdFile(args.inpcrd)


    system = prmtop.createSystem(
        nonbondedMethod=NoCutoff,
        constraints=HBonds,
        implicitSolvent=OBC2,
        implicitSolventSaltConc=args.salt_molar * unit.molar
    )

    dt = args.dt_fs * unit.femtoseconds
    integrator = mm.LangevinIntegrator(
        args.t_start * unit.kelvin,
        args.friction / unit.picosecond,
        dt
    )
    integrator.setRandomNumberSeed(args.seed)

    platform, platform_props = pick_platform(args.platform)
    simulation = Simulation(prmtop.topology, system, integrator, platform, platform_props)
    simulation.context.setPositions(inpcrd.positions)


    t0 = time.time()

    simulation.context.setVelocitiesToTemperature(args.t_start * unit.kelvin, args.seed)
    simulation.minimizeEnergy(
        tolerance=args.min_tol_kj_per_mol_nm * unit.kilojoule_per_mole / unit.nanometer,
        maxIterations=args.min_iter
    )

    total_heat_steps = ps_to_steps(args.heat_ps, args.dt_fs)
    deltaT = max(0.0, args.t_target - args.t_start)

    if deltaT == 0.0 or total_heat_steps == 0:
        integrator.setTemperature(args.t_target * unit.kelvin)
    else:
        n_chunks = int(round(deltaT))  # ~1 K per chunk
        n_chunks = max(1, n_chunks)

        steps_per_chunk = max(1, total_heat_steps // n_chunks)

        simulation.context.setVelocitiesToTemperature(args.t_start * unit.kelvin, args.seed)

        steps_done = 0
        for i in range(n_chunks):
            frac = (i + 1) / float(n_chunks)
            T = args.t_start + frac * (args.t_target - args.t_start)
            integrator.setTemperature(T * unit.kelvin)

            simulation.step(steps_per_chunk)
            steps_done += steps_per_chunk

        remaining = total_heat_steps - steps_done
        if remaining > 0:
            integrator.setTemperature(args.t_target * unit.kelvin)
            simulation.step(remaining)

    state = simulation.context.getState(getPositions=True)
    with open(args.out_pdb, "w") as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)

    integrator.setTemperature(args.t_target * unit.kelvin)

    prod_steps = ns_to_steps(args.prod_ns, args.dt_fs)

    frame_interval_steps = ns_to_steps(args.frame_interval_ns, args.dt_fs)
    frame_interval_steps = max(1, frame_interval_steps)

    simulation.reporters.append(DCDReporter(args.out_dcd, frame_interval_steps))

    simulation.reporters.append(StateDataReporter(
        sys.stdout, 5000,
        step=True, time=True,
        potentialEnergy=True, temperature=True,
        progress=True, remainingTime=True, speed=True,
        totalSteps=prod_steps,
        separator="\t"
    ))

    simulation.step(prod_steps)

    total_sec = time.time() - t0
    with open(args.time_txt, "w") as f:
        f.write(f"Total wall time (s): {total_sec:.3f}\n")
        f.write(f"Total wall time (min): {total_sec/60.0:.3f}\n")
        f.write(f"Total wall time (h): {total_sec/3600.0:.3f}\n")

    print(f"\nDone.\nWrote: {args.out_pdb}, {args.out_dcd}, {args.time_txt}\n")

if __name__ == "__main__":
    main()
