#=
Kramers-Wannier Duality Defect in 1+1D Ising Model
===================================================

This code implements the Kramers-Wannier (KW) duality defect using ITensor.

The KW duality maps:
    σᶻᵢ σᶻᵢ₊₁ → τˣᵢ     (domain wall to spin flip)
    σˣᵢ → τᶻᵢ₋₁ τᶻᵢ    (spin flip to domain wall)

A duality defect is inserted at site d, connecting original spins (i < d)
to dual spins (i ≥ d) through a topological interface.

At criticality (J = h), the defect has quantum dimension √2.

Author: Claude
Date: 2025-11-23
=#

using ITensors
using Printf

"""
    build_ising_with_kw_defect(N, J, h, defect_pos)

Build the Ising Hamiltonian with a Kramers-Wannier duality defect.

The Hamiltonian is:
- Left of defect (i < d): Standard Ising terms
- At defect: Special coupling connecting original and dual descriptions
- Right of defect (i ≥ d): Dual Ising terms

Parameters:
- N: System size
- J: Ising coupling strength
- h: Transverse field strength
- defect_pos: Position of the duality defect

At the critical point J = h, self-duality is manifest.
"""
function build_ising_with_kw_defect(N::Int, J::Real, h::Real, defect_pos::Int)
    sites = siteinds("S=1/2", N; conserve_qns=false)

    # Build MPO manually
    os = OpSum()

    # Region I: Standard Ising (sites 1 to defect_pos-1)
    for i in 1:(defect_pos-2)
        # ZZ coupling
        os += -J, "Sz", i, "Sz", i+1
        # Transverse field
        os += -h, "Sx", i
    end

    # Last site before defect
    if defect_pos > 1
        os += -h, "Sx", defect_pos-1
    end

    # Defect term: This is the key!
    # The KW defect couples σᶻ_{d-1} to the dual description
    # At the defect, we have: σᶻ_{d-1} σᶻ_d = τˣ_{d-1}
    # But we need to represent this in the original basis

    # Defect coupling: Modified ZZ term that acts as the interface
    # The defect acts as: σᶻ_{d-1} connects to σˣ_d (which is τᶻ in dual)
    os += -J, "Sz", defect_pos-1, "Sz", defect_pos

    # Defect transverse field with √2 factor (from fusion rules)
    # This encodes the quantum dimension of the defect
    os += -h * sqrt(2), "Sx", defect_pos

    # Region II: Dual Ising (sites defect_pos+1 to N)
    # In the dual picture: τᶻᵢ τᶻᵢ₊₁ terms become σˣᵢ
    # and τˣᵢ terms become σᶻᵢ σᶻᵢ₊₁

    # At the critical point, we keep the same form (self-dual)
    for i in (defect_pos+1):(N-1)
        # In dual region, exchange role of J and h for proper duality
        os += -h, "Sz", i, "Sz", i+1
        os += -J, "Sx", i
    end

    # Last site
    if defect_pos < N
        os += -J, "Sx", N
    end

    H = MPO(os, sites)
    return H, sites
end


"""
    build_critical_ising_with_defect(N, defect_pos)

Build the critical Ising model (J = h = 1) with KW duality defect.
At criticality, the model is self-dual under KW transformation.
"""
function build_critical_ising_with_defect(N::Int, defect_pos::Int)
    return build_ising_with_kw_defect(N, 1.0, 1.0, defect_pos)
end


"""
    build_defect_hamiltonian_explicit(N, defect_pos)

Build Hamiltonian with explicit defect operator insertion.

This version explicitly constructs the duality defect as a topological
operator that implements the KW transformation at a point.
"""
function build_defect_hamiltonian_explicit(N::Int, defect_pos::Int)
    sites = siteinds("S=1/2", N; conserve_qns=false)

    os = OpSum()

    # Standard Ising everywhere
    for i in 1:(N-1)
        os += -1.0, "Sz", i, "Sz", i+1
        os += -1.0, "Sx", i
    end
    os += -1.0, "Sx", N

    # Now add the defect as a twist
    # The KW defect at bond (d-1, d) is implemented by:
    # D = exp(iπ/4 * σˣ_d)
    # This creates a branch cut for the duality transformation

    # In terms of Hamiltonian modification, the defect changes
    # the coupling at the interface
    # We modify the ZZ term at the defect position

    # Remove standard term and add defect term
    os += 1.0, "Sz", defect_pos-1, "Sz", defect_pos  # Cancel original

    # Add defect coupling: σᶻ_{d-1} σˣ_d (mixes Z and X)
    # This is the hallmark of the duality defect
    os += -1.0, "Sz", defect_pos-1, "Sx", defect_pos

    H = MPO(os, sites)
    return H, sites
end


"""
    run_dmrg_with_defect(H, sites; chi_max=100, nsweeps=10)

Run DMRG to find ground state of Hamiltonian with defect.
"""
function run_dmrg_with_defect(H, sites; chi_max::Int=100, nsweeps::Int=10)
    N = length(sites)

    # Initialize with random product state
    psi0 = randomMPS(sites, linkdims=10)

    # DMRG parameters
    sweeps = Sweeps(nsweeps)
    setmaxdim!(sweeps, 10, 20, 50, chi_max)
    setcutoff!(sweeps, 1e-10)
    setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)

    # Run DMRG
    energy, psi = dmrg(H, psi0, sweeps; outputlevel=0)

    return psi, energy
end


"""
    calculate_entanglement_entropy(psi)

Calculate entanglement entropy at each bond.
"""
function calculate_entanglement_entropy(psi)
    N = length(psi)
    S = zeros(N-1)

    for b in 1:(N-1)
        orthogonalize!(psi, b)
        _, S_vals, _ = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))

        # von Neumann entropy
        p = [s^2 for s in storage(S_vals)]
        p = p[p .> 1e-15]
        S[b] = -sum(p .* log.(p))
    end

    return S
end


"""
    calculate_defect_entropy(psi, defect_pos)

Calculate the entanglement entropy at the defect position.

For the KW defect at criticality, the entropy should show a
characteristic contribution from the defect quantum dimension.
"""
function calculate_defect_entropy(psi, defect_pos::Int)
    N = length(psi)

    if defect_pos < 1 || defect_pos >= N
        error("Defect position must be between 1 and N-1")
    end

    orthogonalize!(psi, defect_pos)
    _, S_vals, _ = svd(psi[defect_pos],
                      (linkind(psi, defect_pos-1), siteind(psi, defect_pos)))

    p = [s^2 for s in storage(S_vals)]
    p = p[p .> 1e-15]
    S_defect = -sum(p .* log.(p))

    return S_defect
end


"""
    measure_correlations(psi, op1, op2, i, j)

Measure two-point correlation function ⟨op1_i op2_j⟩.
"""
function measure_correlations(psi, op1::String, op2::String, i::Int, j::Int)
    if i > j
        i, j = j, i
    end

    orthogonalize!(psi, i)

    # Build contraction
    C = psi[i] * op(op1, siteind(psi, i))
    for k in (i+1):(j-1)
        C *= psi[k]
    end
    C *= psi[j] * op(op2, siteind(psi, j))

    # Contract with conjugate
    for k in j:-1:i
        C *= dag(prime(psi[k], "Link"))
    end

    return scalar(C)
end


"""
    calculate_order_parameters(psi, defect_pos)

Calculate order parameters on both sides of the defect.
"""
function calculate_order_parameters(psi, defect_pos::Int)
    N = length(psi)

    # Magnetization profile
    Sz = zeros(N)
    Sx = zeros(N)

    for i in 1:N
        orthogonalize!(psi, i)
        Sz[i] = real(scalar(psi[i] * op("Sz", siteind(psi, i)) *
                           dag(prime(psi[i], "Site"))))
        Sx[i] = real(scalar(psi[i] * op("Sx", siteind(psi, i)) *
                           dag(prime(psi[i], "Site"))))
    end

    return Sz, Sx
end


"""
    calculate_defect_g_factor(S_with_defect, S_without_defect, defect_pos)

Calculate the g-factor (ground state degeneracy) of the defect.

g = exp(S_defect - S_bulk)

For the KW defect: g = √2
"""
function calculate_defect_g_factor(S_with, S_without, defect_pos)
    ΔS = S_with[defect_pos] - S_without[defect_pos]
    g = exp(ΔS)
    return g
end


"""
    main()

Main function demonstrating the KW duality defect calculation.
"""
function main()
    println("="^60)
    println("Kramers-Wannier Duality Defect in 1+1D Ising Model")
    println("="^60)
    println()

    # Parameters
    N = 40              # System size
    defect_pos = N ÷ 2  # Defect in the middle
    chi_max = 100       # Bond dimension

    println("Parameters:")
    println("  System size N = $N")
    println("  Defect position = $defect_pos")
    println("  Bond dimension χ = $chi_max")
    println()

    # ========================================
    # 1. Critical Ising without defect (reference)
    # ========================================
    println("-"^60)
    println("1. Reference: Critical Ising without defect")
    println("-"^60)

    sites_ref = siteinds("S=1/2", N; conserve_qns=false)
    os_ref = OpSum()
    for i in 1:(N-1)
        os_ref += -1.0, "Sz", i, "Sz", i+1
        os_ref += -1.0, "Sx", i
    end
    os_ref += -1.0, "Sx", N
    H_ref = MPO(os_ref, sites_ref)

    psi_ref, E_ref = run_dmrg_with_defect(H_ref, sites_ref;
                                          chi_max=chi_max, nsweeps=15)
    S_ref = calculate_entanglement_entropy(psi_ref)

    @printf("  Ground state energy: E = %.10f\n", E_ref)
    @printf("  Energy per site: E/N = %.10f\n", E_ref/N)
    @printf("  Entropy at middle bond: S = %.6f\n", S_ref[N÷2])
    println()

    # ========================================
    # 2. Critical Ising with KW duality defect
    # ========================================
    println("-"^60)
    println("2. Critical Ising with KW duality defect")
    println("-"^60)

    H_defect, sites_defect = build_defect_hamiltonian_explicit(N, defect_pos)

    psi_defect, E_defect = run_dmrg_with_defect(H_defect, sites_defect;
                                                 chi_max=chi_max, nsweeps=15)
    S_defect = calculate_entanglement_entropy(psi_defect)

    @printf("  Ground state energy: E = %.10f\n", E_defect)
    @printf("  Energy per site: E/N = %.10f\n", E_defect/N)
    @printf("  Entropy at defect: S = %.6f\n", S_defect[defect_pos])
    println()

    # ========================================
    # 3. Defect properties
    # ========================================
    println("-"^60)
    println("3. Defect Properties")
    println("-"^60)

    # Calculate g-factor
    g_factor = calculate_defect_g_factor(S_defect, S_ref, defect_pos)
    @printf("  Defect g-factor: g = %.6f\n", g_factor)
    @printf("  Expected (KW): g = √2 = %.6f\n", sqrt(2))
    @printf("  Deviation: |g - √2| = %.2e\n", abs(g_factor - sqrt(2)))
    println()

    # Entropy difference at defect
    ΔS = S_defect[defect_pos] - S_ref[defect_pos]
    @printf("  Entropy contribution: ΔS = %.6f\n", ΔS)
    @printf("  Expected: ln(√2) = %.6f\n", log(sqrt(2)))
    println()

    # ========================================
    # 4. Entanglement entropy profile
    # ========================================
    println("-"^60)
    println("4. Entanglement Entropy Profile")
    println("-"^60)

    println("  Bond    S_ref     S_defect   ΔS")
    println("  " * "-"^40)

    bonds_to_show = [1, N÷4, defect_pos-1, defect_pos, defect_pos+1, 3N÷4, N-1]
    for b in bonds_to_show
        if b >= 1 && b <= N-1
            ΔS_b = S_defect[b] - S_ref[b]
            @printf("  %3d    %.4f    %.4f    %+.4f\n",
                    b, S_ref[b], S_defect[b], ΔS_b)
        end
    end
    println()

    # ========================================
    # 5. Order parameters
    # ========================================
    println("-"^60)
    println("5. Order Parameters")
    println("-"^60)

    Sz_ref, Sx_ref = calculate_order_parameters(psi_ref, defect_pos)
    Sz_def, Sx_def = calculate_order_parameters(psi_defect, defect_pos)

    println("  Magnetization ⟨Sᶻ⟩ profile:")
    @printf("    Left of defect (avg):  %.6f\n",
            mean(Sz_def[1:defect_pos-1]))
    @printf("    At defect:             %.6f\n", Sz_def[defect_pos])
    @printf("    Right of defect (avg): %.6f\n",
            mean(Sz_def[defect_pos+1:N]))
    println()

    # ========================================
    # Summary
    # ========================================
    println("="^60)
    println("Summary")
    println("="^60)
    println()
    println("The Kramers-Wannier duality defect has been implemented.")
    println("Key results:")
    @printf("  • Defect quantum dimension: g ≈ %.3f (expected √2 ≈ 1.414)\n", g_factor)
    @printf("  • Topological entropy contribution: ΔS ≈ %.3f\n",
            S_defect[defect_pos] - S_ref[defect_pos])
    println("  • The defect connects original and dual Ising descriptions")
    println()
    println("="^60)

    return psi_defect, E_defect, S_defect, g_factor
end


# Helper function
function mean(x)
    return sum(x) / length(x)
end


# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
