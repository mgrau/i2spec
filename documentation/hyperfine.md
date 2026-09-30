# Hyperfine structure

## Angular momenta

Both iodine nuclei have spin, i = 5/2 for ¹²⁷I and i = 7/2 for ¹²⁹I, and an electric quadrupole
moment. The nuclear spins couple to the total nuclear spin **I** = **I**₁ + **I**₂, and **I** couples
to the rotational angular momentum **J** to give the total angular momentum **F** = **J** + **I**.
The hyperfine Hamiltonian is set up in the coupled basis |J, (i₁ i₂) I; F⟩ and diagonalised
separately for each F. For ¹²⁷I₂ each rovibrational level has 15 or 21 hyperfine sublevels, depending
on the parity of J (see below).

## Effective hyperfine Hamiltonian

The effective Hamiltonian of Broyer, Vigué & Lehmann (1978) is

$$
H_\text{hfs} = eQq\,H_\text{EQ} + C\,H_\text{SR} + d\,H_\text{TSS} + \delta\,H_\text{SSS}, \tag{1}
$$

with the electric quadrupole interaction (eQq), the nuclear spin–rotation interaction (C), and the
tensor (d) and scalar (δ) nuclear spin–spin interactions. Within one value of J the operators are

$$
H_\text{EQ} = -\sum_{n=1,2}\frac{3(\mathbf I_n\cdot\mathbf J)^2 + \tfrac32(\mathbf I_n\cdot\mathbf J) - \mathbf I_n^2\mathbf J^2}{2i_n(2i_n-1)(2J-1)(2J+3)},\qquad H_\text{SR} = \mathbf I_1\cdot\mathbf J + \frac{C_2}{C_1}\,\mathbf I_2\cdot\mathbf J, \tag{2}
$$

$$
H_\text{TSS} = \frac{3(\mathbf I_1\cdot\mathbf J)(\mathbf I_2\cdot\mathbf J) + 3(\mathbf I_2\cdot\mathbf J)(\mathbf I_1\cdot\mathbf J) - 2(\mathbf I_1\cdot\mathbf I_2)\,J(J+1)}{(2J-1)(2J+3)},\qquad H_\text{SSS} = \mathbf I_1\cdot\mathbf I_2 . \tag{3}
$$

The quadrupole and tensor spin–spin operators are second-rank tensors with respect to the molecular
axis. They therefore also have matrix elements between J and J ± 2 (and between I and I ± 2). These
matrix elements shift the components by up to a few tens of kHz at high J and are included; the
energies of the neighbouring rotational levels are taken from the level model. All matrix elements are
evaluated with Wigner 3j, 6j and 9j symbols. For ¹²⁷I¹²⁹I the two nuclei have different eQq and C,
which gives additional matrix elements with ΔI = ±1.

## Nuclear-spin statistics

In ¹²⁷I₂ and ¹²⁹I₂ the two nuclei are identical fermions. The total wavefunction must be antisymmetric
under their exchange, so in X (0g⁺) only values of I with the same parity as J occur, and in
B (0u⁺) only values of the opposite parity. For ¹²⁷I₂ (I = 0 … 5) the nuclear-spin statistical weight
is therefore 15 for even J″ (I = 0, 2, 4) and 21 for odd J″ (I = 1, 3, 5); a line with even J″ has 15
main components and a line with odd J″ has 21.

<figure class="fig" markdown="span">
<div class="svg" data-svg="figures/hyperfine_levels.svg" role="img" aria-label="Hyperfine energy levels of the upper and lower levels of R(56) 32-0 and its 15 main transitions"></div>
<figcaption markdown="span">**Figure 1.** Hyperfine sublevels of the lower level X(v″ = 0, J″ = 56) and the upper level B(v′ = 32, J′ = 57) of R(56) 32–0, and the 15 main transitions a1–a15. Energies are relative to the centre of each group. The upper group is drawn with a three times larger vertical scale; the scale bars give 100 MHz (upper) and 500 MHz (lower), and the separation of the two groups is not to scale. The splitting of both levels is dominated by the quadrupole interaction; eQq is −2453 MHz for the lower level and −545 MHz for the upper level, so the lower level is split about 4.5 times more widely (1.1 GHz against 0.25 GHz).</figcaption>
</figure>

## Parameters

The parameters eQq, C, d and δ of each level are taken from interpolation formulae in the
vibrational energy: those of Bodermann, Knöckel & Tiemann (2002) for ¹²⁷I₂ up to v′ = 43, and those of
Salumbides *et al.* (2006) above v′ = 43 and for the other isotopologues, scaled by the ratios of the
nuclear moments. For ¹²⁷I₂ the formulae are supplemented by a table of B-state corrections fitted to
the measured hyperfine spectra at each v′ for which measurements exist.

## Components and their intensities

The relative intensity of a component is the square of the electric-dipole matrix element between the
upper and lower hyperfine eigenvectors, normalised to unity over the line. The dipole operator acts on
the rotational coordinates only, so ΔF = 0, ±1 and ΔI = 0 in the uncoupled limit. The main components
have ΔF = ΔJ and carry almost all of the intensity at high J. They are labelled a1, a2, … in order of
increasing frequency, as in the BIPM tables.

<figure class="fig" markdown="span">
<div class="svg" data-svg="figures/hyperfine.svg" role="img" aria-label="Hyperfine components of R(56) 32-0 and P(53) 32-0 under the Doppler profile"></div>
<figcaption markdown="span">**Figure 2.** Main hyperfine components of R(56) 32–0 (even J″, 15 components) and P(53) 32–0 (odd J″, 21 components), with the Doppler profile at 20 °C. The full width of the pattern, about 1 GHz, is comparable to the Doppler width, so the components are resolved only by sub-Doppler methods.</figcaption>
</figure>

!!! note "Accuracy of the hyperfine offsets"
    The component offsets relative to the hyperfine-free position are accurate to a few tens of kHz
    for ¹²⁷I₂ with v′ ≤ 43, and to about 1 MHz for ¹²⁹I₂ and ¹²⁷I¹²⁹I.
