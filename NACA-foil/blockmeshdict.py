#!/usr/bin/env python

from __future__ import division, print_function
import argparse
import numpy as np
from numpy import linspace, zeros, ones, sin, cos, arctan, pi
import os

def gen_blockmeshdict(foil="0012", alpha_deg=4):
    """
    Write a `blockMeshDict` for a NACA foil at specified angle of attack.
    """
    # Foil geometry
    c = 1.0              # Geometric chord length
    alpha = np.deg2rad(alpha_deg)  # Angle of attack (in radians)
    NACA = [int(d) for d in foil]  # NACA 4-digit designation


    # Mesh dimensions
    scale = 1            # Scaling factor
    H = 8                # *Half* height of channel
    W = 0.5              # *Half* depth of foil (z-direction)
    D = 16               # Length of downstream section

    # Mesh resolution parameters
    Ni = 400             # Number of interpolation points along the foil

    # Nx = 200           # Number of mesh cells along the foil
    Nleading = 40        # Number of mesh cells along the leading foil
    Ntrailing = 20       # Number of mesh cells along the trailing foil
    
    ND = 150             # Number of cells in the downstream direction
    NT = 100             # Number of cells the transverse direction
    NW = 1               # Number of cells in the z-direction (along the foil axis)

    # Expansion rates
    ExpTransverse = 100 # Expansion rate in transverse direction
    ExpDownstream = 1   # Expansion rate in the downstream direction
    ExpLeading = 1      # Expansion rate in the leading foil
    ExpTrailing = 1     # Expansion rate in the trailing foil

    # ------------------------- END OF MESH PARAMETER REGION --------------------- #


    # Create a vector with x-coordinates, camber and thickness
    beta = linspace(0, pi, Ni)
    x = c*(0.5*(1 - cos(beta)))
    y_c = zeros(len(x))
    y_t = zeros(len(x))
    theta = zeros(len(x))


    # Values of m, p and t
    m = NACA[0]/100
    p = NACA[1]/10
    t = (NACA[2]*10 + NACA[3])/100


    # Calculate thickness
    # The upper expression will give the airfoil a finite thickness at the trailing
    # edge, witch might cause trouble. The lower expression is corrected to give
    # zero thickness at the trailing edge, but the foil is strictly speaking no
    # longer a proper NACA airfoil.
    #
    # See http://turbmodels.larc.nasa.gov/naca4412sep_val.html
    #     http://en.wikipedia.org/wiki/NACA_airfoil

    #y_t = (t*c/0.2) * (0.2969*(x/c)**0.5 - 0.1260*(x/c) - 0.3516*(x/c)**2 + 0.2843*(x/c)**3 - 0.1015*(x/c)**4)
    y_t = (t*c/0.2)*(0.2969*(x/c)**0.5 - 0.1260*(x/c) - 0.3516*(x/c)**2 \
           + 0.2843*(x/c)**3 - 0.1036*(x/c)**4)

    if p > 0:
        # Calculate camber
        y_c += (m*x/p**2)*(2*p - x/c)*(x < p*c)
        y_c += (m*(c-x)/(1 - p)**2)*(1 + x/c - 2*p)*(x >= p*c)
        # Calculate theta-value
        theta += arctan((m/p**2) * (2*p - 2*x/c))*(x < p*c)
        theta += arctan((m/(1 - p)**2) * (-2*x/c + 2*p))*(x >= p*c)


    # Calculate coordinates of upper surface
    Xu = x - y_t*sin(theta) - scale/4.0
    Yu = y_c + y_t*cos(theta)


    # Calculate coordinates of lower surface
    Xl = x + y_t*sin(theta) - scale/4.0
    Yl = y_c - y_t*cos(theta)


    # Rotate foil to reach specified angle of attack
    upper = np.matrix([[cos(alpha), sin(alpha)],
                       [-sin(alpha), cos(alpha)]])*np.vstack((Xu, Yu))
    lower = np.matrix([[cos(alpha), sin(alpha)],
                      [-sin(alpha), cos(alpha)]])*np.vstack((Xl, Yl))

    Xu = upper[0, :].conj().transpose()
    Yu = upper[1, :].conj().transpose()
    Xl = lower[0, :].conj().transpose()
    Yl = lower[1, :].conj().transpose()

    if p > 0:
        # Find index i of max. camber
        C_max_idx = np.where(y_c == max(y_c))[0][0]
    else:
        # Otherwise use location of max. thickness
        C_max_idx = np.where(y_t == max(y_t))[0][0]


    # Move point of mesh "nose"
    NoseX = (-H + Xu[C_max_idx])*cos(alpha)
    NoseY = -(-H + Xu[C_max_idx])*sin(alpha)


    # Calculate the location of the vertices on the positive y-axis and put them in a matrix
    vertices = zeros((12, 3))

    vertices[0, :] = [NoseX[0], NoseY[0], W]
    vertices[1, :] = [Xu[C_max_idx], H, W]
    vertices[2, :] = [Xu[-1], H, W]
    vertices[3, :] = [D, H, W]
    vertices[4, :] = [Xu[0], Yu[0], W]
    vertices[5, :] = [Xu[C_max_idx], Yu[C_max_idx], W]
    vertices[6, :] = [Xl[C_max_idx], Yl[C_max_idx], W]
    vertices[7, :] = [Xu[-1], Yu[-1], W]
    vertices[8, :] = [D, Yu[-1], W]
    vertices[9, :] = [Xl[C_max_idx], -H, W]
    vertices[10, :] = [Xu[-1], -H, W]
    vertices[11, :] = [D, -H, W]

    # Create vertices for other side (negative z-axis)
    vertices2 = vertices.copy()
    vertices2[:, 2] *= -1
    vertices = np.vstack((vertices, vertices2))


    # Edge 4-5 and 16-17
    pts1 = np.concatenate([Xu[1:C_max_idx], 
                           Yu[1:C_max_idx], 
                           W*ones(np.shape(Xu[1:C_max_idx]))], axis=1)

    pts5 = np.concatenate([pts1[:, 0], pts1[:, 1], -pts1[:, 2]], axis=1)

    # Edge 5-7 and 17-19
    pts2 = np.concatenate([Xu[C_max_idx + 1:Ni - 1],
                           Yu[C_max_idx + 1:Ni - 1],
                           W*ones(np.shape(Xu[C_max_idx + 1:Ni - 1]))], axis=1)

    pts6 = np.concatenate([pts2[:, 0], pts2[:, 1], -pts2[:, 2]], axis=1)

    # Edge 4-6 and 16-18
    pts3 = np.concatenate([Xl[1:C_max_idx], 
                           Yl[1:C_max_idx],
                           W*ones(np.shape(Xl[1:C_max_idx]))], axis=1)

    pts7 = np.concatenate([pts3[:, 0], pts3[:, 1], -pts3[:, 2]], axis=1)

    # Edge 6-7 and 18-19
    pts4 = np.concatenate([Xl[C_max_idx + 1:Ni - 1],
                           Yl[C_max_idx + 1:Ni - 1],
                           W*ones(np.shape(Xl[C_max_idx + 1:Ni - 1]))], axis=1)
    
    pts8 = np.concatenate([pts4[:, 0], pts4[:, 1], -pts4[:, 2]], axis=1)

    # Edge 0-1 and 12-13
    pts9 = np.array([-H*cos(pi/4) + Xu[C_max_idx, 0], H*sin(pi/4), W])
    pts11 = np.array([pts9[0], pts9[1], -pts9[2]])

    # Edge 0-9 and 12-21
    pts10 = np.array([-H*cos(pi/4) + Xu[C_max_idx, 0], -H*sin(pi/4), W])
    pts12 = np.array([pts10[0], pts10[1], -pts10[2]])

    # Calculate number of mesh points along 4-5 and 4-6
    # Nleading = (C_max_idx/Ni)*Nx
    # Nleading = int((x[C_max_idx]/c)*Nx)
    
    # # Calculate number of mesh points along 5-7 and 6-7
    # Ntrailing = (Nx - Nleading)-130

    # Open file
    f = open("blockMeshDict", "w")


    # Write file
    f.write("/*--------------------------------*- C++ -*----------------------------------*\\ \n")
    f.write("| =========                 |                                                 | \n")
    f.write("| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n")
    f.write("|  \\\\    /   O peration     | Version:  3.0.x                                 | \n")
    f.write("|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      | \n")
    f.write("|    \\\\/     M anipulation  |                                                 | \n")
    f.write("\\*---------------------------------------------------------------------------*/ \n")
    f.write("FoamFile                                                                        \n")
    f.write("{                                                                               \n")
    f.write("    version     2.0;                                                            \n")
    f.write("    format      ascii;                                                          \n")
    f.write("    class       dictionary;                                                     \n")
    f.write("    object      blockMeshDict;                                                  \n")
    f.write("}                                                                               \n")
    f.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n")
    f.write("\n")
    f.write("convertToMeters %f; \n" % scale)
    f.write("\n")
    f.write("vertices \n")
    f.write("( \n")
    for vertex in vertices:
        f.write("    (%f %f %f)\n" % tuple(vertex))
    f.write("); \n")
    f.write("\n")
    f.write("blocks \n")
    f.write("( \n")
    f.write("    hex (16 17 13 12 4 5 1 0)     (%i %i %i) simpleGrading (%f %f 1) \n" % (Nleading, NT, NW, ExpLeading, ExpTransverse))
    f.write("    hex (17 19 14 13 5 7 2 1)     (%i %i %i) simpleGrading (%f %f 1) \n" % (Ntrailing, NT, NW, ExpTrailing, ExpTransverse))
    f.write("    hex (19 20 15 14 7 8 3 2)     (%i %i %i) simpleGrading (%f %f 1) \n" % (ND, NT, NW, ExpDownstream, ExpTransverse))
    f.write("    hex (4 6 9 0 16 18 21 12)     (%i %i %i) simpleGrading (%f %f 1) \n" % (Nleading, NT, NW, ExpLeading, ExpTransverse))
    f.write("    hex (6 7 10 9 18 19 22 21)    (%i %i %i) simpleGrading (%f %f 1) \n" % (Ntrailing, NT, NW, ExpTrailing, ExpTransverse))
    f.write("    hex (7 8 11 10 19 20 23 22)   (%i %i %i) simpleGrading (%f %f 1) \n" % (ND, NT, NW, ExpDownstream, ExpTransverse))

    f.write("); \n")
    f.write("\n")
    f.write("edges \n")
    f.write("( \n")

    f.write("    spline 4 5 \n")
    f.write("        ( \n")
    for pt in np.array(pts1):
        f.write("            (%f %f %f) \n" % tuple(pt))
    f.write("        ) \n")

    f.write("    spline 5 7 \n")
    f.write("        ( \n")
    for pt in np.array(pts2):
        f.write("            (%f %f %f)\n" % tuple(pt))
    f.write("        ) \n")

    f.write("    spline 4 6 \n")
    f.write("        ( \n")
    for pt in np.array(pts3):
        f.write("            (%f %f %f)\n" % tuple(pt))
    f.write("        ) \n")

    f.write("    spline 6 7 \n")
    f.write("        ( \n")
    for pt in np.array(pts4):
        f.write("            (%f %f %f)\n" % tuple(pt))
    f.write("        ) \n")

    f.write("    spline 16 17 \n")
    f.write("        ( \n")
    for pt in np.array(pts5):
        f.write("            (%f %f %f)\n" % tuple(pt))
    f.write("        ) \n")

    f.write("    spline 17 19 \n")
    f.write("        ( \n")
    for pt in np.array(pts6):
        f.write("            (%f %f %f)\n" % tuple(pt))
    f.write("        ) \n")

    f.write("    spline 16 18 \n")
    f.write("        ( \n")
    for pt in np.array(pts7):
        f.write("            (%f %f %f)\n" % tuple(pt))
    f.write("        ) \n")

    f.write("    spline 18 19 \n")
    f.write("        ( \n")
    for pt in np.array(pts8):
        f.write("            (%f %f %f)\n" % tuple(pt))
    f.write("        ) \n")

    f.write("    arc 0 1 (%f %f %f) \n" % tuple(pts9))
    f.write("    arc 0 9 (%f %f %f) \n" % tuple(pts10))
    f.write("    arc 12 13 (%f %f %f) \n" % tuple(pts11))
    f.write("    arc 12 21 (%f %f %f) \n" % tuple(pts12))

    f.write("); \n")
    f.write("\n")
    f.write("boundary \n")
    f.write("( \n")

    f.write("    inlet \n")
    f.write("    { \n")
    f.write("        type patch; \n")
    f.write("        faces \n")
    f.write("        ( \n")
    f.write("            (1 0 12 13) \n")
    f.write("            (0 9 21 12) \n")
    f.write("        ); \n")
    f.write("    } \n")
    f.write("\n")

    f.write("    outlet \n")
    f.write("    { \n")
    f.write("        type patch; \n")
    f.write("        faces \n")
    f.write("        ( \n")
    f.write("            (11 8 20 23) \n")
    f.write("            (8 3 15 20) \n")
    f.write("        ); \n")
    f.write("    } \n")
    f.write("\n")

    f.write("    topAndBottom \n")
    f.write("    { \n")
    f.write("        type patch; \n")
    f.write("        faces \n")
    f.write("        ( \n")
    f.write("            (3 2 14 15) \n")
    f.write("            (2 1 13 14) \n")
    f.write("            (9 10 22 21) \n")
    f.write("            (10 11 23 22) \n")
    f.write("        ); \n")
    f.write("    } \n")
    f.write("\n")

    f.write("    airfoil \n")
    f.write("    { \n")
    f.write("        type wall; \n")
    f.write("        faces \n")
    f.write("        ( \n")
    f.write("            (5 4 16 17) \n")
    f.write("            (7 5 17 19) \n")
    f.write("            (4 6 18 16) \n")
    f.write("            (6 7 19 18) \n")
    f.write("        ); \n")
    f.write("    } \n")
    f.write("); \n")
    f.write(" \n")
    f.write("mergePatchPairs \n")
    f.write("( \n")
    f.write("); \n")
    f.write(" \n")
    f.write("// ************************************************************************* // \n")

    # Close file
    f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plotting results")
    parser.add_argument("foil", help="NACA foil digits")
    parser.add_argument("alpha_deg", help="Angle of attack (deg)")
    args = parser.parse_args()

    print("Generating blockMeshDict for a NACA {} at {} "
          "degrees angle of attack".format(args.foil, args.alpha_deg))

    gen_blockmeshdict(args.foil, float(args.alpha_deg))
