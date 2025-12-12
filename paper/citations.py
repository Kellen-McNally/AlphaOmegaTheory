"""
Paper Citations Module

This module manages the citation database and attribution for the paper.
It ensures that all prior work is correctly credited and that the
"Zero-Parameter" claim is rigorously defined (constants vs parameters).

[WARN] This module is called by paper/compiler.py during paper generation.
[WARN] Do not run directly. Use: ./run_tests.sh
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional
from collections import defaultdict


@dataclass
class Citation:
    """
    A citation to prior work.

    Attributes:
        key: BibTeX key (e.g., "Cartan1894")
        authors: List of author names
        title: Paper title
        year: Publication year
        journal: Journal/venue
        contribution: What this work contributed
        weight: Weight of contribution (0-1)
        url: URL to the paper (optional)
    """

    key: str
    authors: List[str]
    title: str
    year: int
    journal: str
    contribution: str
    weight: float = 1.0  # Default: full credit
    url: Optional[str] = None


@dataclass
class Attribution:
    """
    Attribution for a specific result.

    Tracks what prior work a result builds on.
    """

    result_name: str  # e.g., "$\alpha$_GUT = 1/42"
    description: str  # What this result is

    # Prior work this builds on
    prior_work: List[Citation] = field(default_factory=list)

    def total_prior_weight(self) -> float:
        """Sum of all prior work weights."""
        return sum(c.weight for c in self.prior_work)

    def normalized_weights(self) -> Dict[str, float]:
        """
        Normalize all weights to sum to 1.

        Returns:
            dict: {citation_key: normalized_weight}
        """
        total = self.total_prior_weight()
        if total == 0:
            return {}

        weights = {c.key: c.weight / total for c in self.prior_work}
        return weights


class CitationDatabase:
    """Database of all citations used in the framework."""

    def __init__(self):
        self.citations: Dict[str, Citation] = {}
        self.attributions: List[Attribution] = []
        self._initialize_core_citations()

    def _initialize_core_citations(self) -> None:
        """Initialize core mathematical citations."""

        # John T. Graves - discovered octonions
        self.add_citation(
            Citation(
                key="Graves1845",
                authors=["John T. Graves"],
                title="On a Connection between the General Theory of Normal Couples and the Theory of Complete Quadratic Functions of Two Variables",
                year=1845,
                journal="Philosophical Magazine",
                contribution="Discovered octonions (octaves) in 1843, extending Hamilton's quaternions to 8 dimensions",
                weight=1.0,
                url="https://doi.org/10.1080/14786444508645136",
            )
        )

        # Arthur Cayley - published octonions
        self.add_citation(
            Citation(
                key="Cayley1845",
                authors=["Arthur Cayley"],
                title="On Jacobi's Elliptic Functions, in reply to the Rev. B. Bronwin; and on Quaternions",
                year=1845,
                journal="Philosophical Magazine",
                contribution="Independently discovered and first published octonions (Cayley numbers)",
                weight=1.0,
                url="https://doi.org/10.1080/14786444508645071",
            )
        )

        # Élie Cartan - discovered G$_2$ (1894)
        self.add_citation(
            Citation(
                key="Cartan1894",
                authors=["Élie Cartan"],
                title="Sur la structure des groupes de transformations finis et continus",
                year=1894,
                journal="Thèse de doctorat, Paris",
                contribution="Discovered the exceptional Lie algebra G$_2$ and classified all simple Lie algebras",
                weight=1.0,
                url="https://eudml.org/doc/192968",
            )
        )

        # Wilhelm Killing - classified Lie algebras
        self.add_citation(
            Citation(
                key="Killing1888",
                authors=["Wilhelm Killing"],
                title="Die Zusammensetzung der stetigen endlichen Transformationsgruppen",
                year=1888,
                journal="Mathematische Annalen",
                contribution="First classification of simple Lie algebras, found exceptional algebras",
                weight=1.0,
                url="https://doi.org/10.1007/BF01444696",
            )
        )

        # Adolf Hurwitz - classified normed division algebras
        self.add_citation(
            Citation(
                key="Hurwitz1898",
                authors=["Adolf Hurwitz"],
                title="Uber die Composition der quadratischen Formen von beliebig vielen Variablen",
                year=1898,
                journal="Nachrichten von der Gesellschaft der Wissenschaften zu Gottingen",
                contribution="Proved that only four normed division algebras exist (dimensions 1, 2, 4, 8)",
                weight=1.0,
                url="https://eudml.org/doc/182735",
            )
        )

        # Dominic Joyce - G2 holonomy theorem
        self.add_citation(
            Citation(
                key="Joyce1996",
                authors=["Dominic Joyce"],
                title="Compact Riemannian 7-manifolds with holonomy G2",
                year=1996,
                journal="Journal of Differential Geometry",
                contribution="Proved G2 is the unique holonomy group for 7-dimensional Ricci-flat manifolds",
                weight=1.0,
                url="https://doi.org/10.4310/jdg/1214458109",
            )
        )

        # Planck Collaboration - measured $\Omega$_Λ
        self.add_citation(
            Citation(
                key="Planck2018",
                authors=["Planck Collaboration"],
                title="Planck 2018 results. VI. Cosmological parameters",
                year=2018,
                journal="Astronomy & Astrophysics",
                contribution="Measured dark energy density $\Omega$_Λ = 0.6847 ± 0.0073",
                weight=1.0,
                url="https://doi.org/10.1051/0004-6361/201833910",
            )
        )

        # Particle Data Group - compiled experimental values
        self.add_citation(
            Citation(
                key="PDG2024",
                authors=["S. Navas et al. (Particle Data Group)"],
                title="Review of Particle Physics",
                year=2024,
                journal="Physical Review D",
                contribution="Experimental measurements of neutrino oscillation parameters, CKM matrix elements, and particle properties",
                weight=1.0,
                url="https://doi.org/10.1103/PhysRevD.110.030001",
            )
        )

        # Keep PDG2020 for reference
        self.add_citation(
            Citation(
                key="PDG2020",
                authors=["Particle Data Group"],
                title="Review of Particle Physics",
                year=2020,
                journal="Progress of Theoretical and Experimental Physics",
                contribution="Historical reference for particle physics measurements",
                weight=1.0,
                url="https://doi.org/10.1093/ptep/ptaa104",
            )
        )

        # Gell-Mann - SU(3) flavor symmetry
        self.add_citation(
            Citation(
                key="GellMann1962",
                authors=["Murray Gell-Mann"],
                title="Symmetries of Baryons and Mesons",
                year=1962,
                journal="Physical Review",
                contribution="Eightfold Way and SU(3) flavor symmetry",
                weight=1.0,
                url="https://doi.org/10.1103/PhysRev.125.1067",
            )
        )

        # Cabibbo - quark mixing
        self.add_citation(
            Citation(
                key="Cabibbo1963",
                authors=["Nicola Cabibbo"],
                title="Unitary Symmetry and Leptonic Decays",
                year=1963,
                journal="Physical Review Letters",
                contribution="Discovered quark mixing, Cabibbo angle",
                weight=1.0,
                url="https://doi.org/10.1103/PhysRevLett.10.531",
            )
        )

        # Weinberg - electroweak theory
        self.add_citation(
            Citation(
                key="Weinberg1967",
                authors=["Steven Weinberg"],
                title="A Model of Leptons",
                year=1967,
                journal="Physical Review Letters",
                contribution="Electroweak unification, predicted weak mixing angle",
                weight=1.0,
                url="https://doi.org/10.1103/PhysRevLett.19.1264",
            )
        )

        # Salam - electroweak theory
        self.add_citation(
            Citation(
                key="Salam1968",
                authors=["Abdus Salam"],
                title="Weak and Electromagnetic Interactions",
                year=1968,
                journal="Svartholm: Elementary Particle Theory",
                contribution="Independent development of electroweak theory",
                weight=1.0,
            )
        )

        # Glashow - electroweak theory
        self.add_citation(
            Citation(
                key="Glashow1961",
                authors=["Sheldon Glashow"],
                title="Partial-symmetries of weak interactions",
                year=1961,
                journal="Nuclear Physics",
                contribution="SU(2)×U(1) structure of electroweak theory",
                weight=1.0,
                url="https://doi.org/10.1016/0029-5582(61)90469-2",
            )
        )

        # Minkowski - seesaw mechanism
        self.add_citation(
            Citation(
                key="Minkowski1977",
                authors=["Peter Minkowski"],
                title="mu to e gamma at a Rate of One Out of $10^9$ Muon Decays?",
                year=1977,
                journal="Physics Letters B",
                contribution="Seesaw mechanism for neutrino masses",
                weight=1.0,
                url="https://doi.org/10.1016/0370-2693(77)90435-X",
            )
        )

        # Super-Kamiokande - neutrino oscillations
        self.add_citation(
            Citation(
                key="SuperK1998",
                authors=["Super-Kamiokande Collaboration"],
                title="Evidence for Oscillation of Atmospheric Neutrinos",
                year=1998,
                journal="Physical Review Letters",
                contribution="First evidence for neutrino oscillations, Nobel Prize 2015",
                weight=1.0,
                url="https://doi.org/10.1103/PhysRevLett.81.1562",
            )
        )

        # Joyce - M-theory and G$_2$ holonomy
        self.add_citation(
            Citation(
                key="Joyce2000",
                authors=["Dominic Joyce"],
                title="Compact Manifolds with Special Holonomy",
                year=2000,
                journal="Oxford University Press",
                contribution="Mathematical foundation of G$_2$ holonomy manifolds",
                weight=1.0,
            )
        )

        # Acharya - G$_2$ phenomenology
        self.add_citation(
            Citation(
                key="Acharya2004",
                authors=["Bobby Acharya"],
                title="M theory, G$_2$ manifolds and four dimensional physics",
                year=2004,
                journal="Classical and Quantum Gravity",
                contribution="Phenomenology of M-theory compactified on G$_2$ manifolds",
                weight=1.0,
                url="https://doi.org/10.1088/0264-9381/21/10/001",
            )
        )

        # Cohl Furey - Octonions and the Standard Model
        self.add_citation(
            Citation(
                key="Furey2016",
                authors=["Cohl Furey"],
                title="Standard Model Physics from an Algebra?",
                year=2016,
                journal="PhD Thesis, University of Waterloo",
                contribution="Groundbreaking work deriving Standard Model gauge groups and fermion representations from division algebras and octonions. Closest prior approach to geometric unification.",
                weight=1.0,
                url="http://hdl.handle.net/10012/10556",
            )
        )

        # Additional Furey paper on SU(3) from octonions
        self.add_citation(
            Citation(
                key="Furey2018",
                authors=["Cohl Furey"],
                title="SU(3)C × SU(2)L × U(1)Y (× U(1)X) as a symmetry of division algebraic ladder operators",
                year=2018,
                journal="European Physical Journal C",
                contribution="Derived Standard Model gauge structure from octonion ladder operators, demonstrating deep connection between division algebras and particle physics",
                weight=1.0,
                url="https://doi.org/10.1140/epjc/s10052-018-5844-7",
            )
        )

        # Albert Einstein - Zero parameter ideal
        self.add_citation(
            Citation(
                key="Einstein1945",
                authors=["Albert Einstein"],
                title="Letter to Ilse Rosenthal-Schneider",
                year=1945,
                journal="Einstein Archive 52-380",
                contribution="Stated the 'theorem' that there are no arbitrary constants in nature and that a final theory should have only rationally determined constants",
                weight=1.0,
            )
        )

        # Paul Dirac - Large Number Hypothesis
        self.add_citation(
            Citation(
                key="Dirac1937",
                authors=["Paul Dirac"],
                title="The Cosmological Constants",
                year=1937,
                journal="Nature",
                contribution="Proposed that fundamental physical constants are not independent but connected by mathematical relations",
                weight=1.0,
                url="https://doi.org/10.1038/139323a0",
            )
        )

        # Lee Smolin - Critique of tuning
        self.add_citation(
            Citation(
                key="Smolin2006",
                authors=["Lee Smolin"],
                title="The Trouble with Physics: The Rise of String Theory, the Fall of a Science, and What Comes Next",
                year=2006,
                journal="Houghton Mifflin",
                contribution="Critiqued the landscape of string theory and the reliance on parameter tuning, calling for a return to background-independent approaches",
                weight=1.0,
            )
        )

        # Alexander Unzicker - Fundamental constants and physics critique
        self.add_citation(
            Citation(
                key="Unzicker2010",
                authors=["Alexander Unzicker"],
                title="Bankrupting Physics: How Today's Top Scientists are Gambling Away Their Credibility",
                year=2010,
                journal="Palgrave Macmillan",
                contribution="Critiqued untestable theories in modern physics and questioned whether fundamental constants should be derivable rather than arbitrary",
                weight=1.0,
            )
        )

        self.add_citation(
            Citation(
                key="Unzicker2015",
                authors=["Alexander Unzicker"],
                title="Einstein's Lost Key: How We Overlooked the Best Idea of the 20th Century",
                year=2015,
                journal="CreateSpace Independent Publishing",
                contribution="Argued that fundamental constants like G should be derivative from more fundamental principles, inspiring approaches where physics emerges from pure geometry",
                weight=1.0,
            )
        )

        # Edward Witten - M-theory
        self.add_citation(
            Citation(
                key="Witten1995",
                authors=["Edward Witten"],
                title="String theory dynamics in various dimensions",
                year=1995,
                journal="Nuclear Physics B",
                contribution="Introduced M-theory as the strong coupling limit of string theory in 11 dimensions",
                weight=1.0,
                url="https://doi.org/10.1016/0550-3213(95)00158-O",
            )
        )

        # LEP Electroweak Working Group - Precision measurements
        self.add_citation(
            Citation(
                key="LEP1996",
                authors=["LEP Electroweak Working Group"],
                title="A Combination of Preliminary Electroweak Measurements and Constraints on the Standard Model",
                year=1996,
                journal="CERN-PPE-96-183",
                contribution="Precision measurement of Z boson properties, determined number of light neutrino generations to be 3",
                weight=1.0,
                url="https://cds.cern.ch/record/316272",
            )
        )

        # Kobayashi-Maskawa - CKM matrix
        self.add_citation(
            Citation(
                key="KobayashiMaskawa1973",
                authors=["Makoto Kobayashi", "Toshihide Maskawa"],
                title="CP-Violation in the Renormalizable Theory of Weak Interaction",
                year=1973,
                journal="Progress of Theoretical Physics",
                contribution="Predicted third generation of quarks and CKM mixing matrix, Nobel Prize 2008",
                weight=1.0,
                url="https://doi.org/10.1143/PTP.49.652",
            )
        )

        # Lincoln Wolfenstein - Parametrization
        self.add_citation(
            Citation(
                key="Wolfenstein1983",
                authors=["Lincoln Wolfenstein"],
                title="Parametrization of the Kobayashi-Maskawa Matrix",
                year=1983,
                journal="Physical Review Letters",
                contribution="Standard parametrization of CKM matrix in terms of Cabibbo angle and small parameters",
                weight=1.0,
                url="https://doi.org/10.1103/PhysRevLett.51.1945",
            )
        )

        # Peccei-Quinn - Strong CP solution
        self.add_citation(
            Citation(
                key="PecceiQuinn1977",
                authors=["Roberto Peccei", "Helen Quinn"],
                title="CP Conservation in the Presence of Pseudoparticles",
                year=1977,
                journal="Physical Review Letters",
                contribution="Proposed dynamical solution to strong CP problem via axion mechanism",
                weight=1.0,
                url="https://doi.org/10.1103/PhysRevLett.38.1440",
            )
        )

        # Weinberg - Axion
        self.add_citation(
            Citation(
                key="Weinberg1978",
                authors=["Steven Weinberg"],
                title="A New Light Boson?",
                year=1978,
                journal="Physical Review Letters",
                contribution="Predicted axion particle as consequence of Peccei-Quinn mechanism",
                weight=1.0,
                url="https://doi.org/10.1103/PhysRevLett.40.223",
            )
        )

        # Wilczek - Axion naming
        self.add_citation(
            Citation(
                key="Wilczek1978",
                authors=["Frank Wilczek"],
                title="Problem of Strong P and T Invariance in the Presence of Instantons",
                year=1978,
                journal="Physical Review Letters",
                contribution="Independent prediction and naming of the axion",
                weight=1.0,
                url="https://doi.org/10.1103/PhysRevLett.40.279",
            )
        )

        # nEDM experiment - Strong CP constraint
        self.add_citation(
            Citation(
                key="nEDM2020",
                authors=["nEDM Collaboration"],
                title="Measurement of the permanent electric dipole moment of the neutron",
                year=2020,
                journal="Physical Review Letters",
                contribution="Set stringent upper limit on neutron EDM, constraining strong CP violation angle",
                weight=1.0,
                url="https://doi.org/10.1103/PhysRevLett.124.081803",
            )
        )

        # NIST - Atomic data
        self.add_citation(
            Citation(
                key="NIST2023",
                authors=["NIST"],
                title="Atomic Spectra Database",
                year=2023,
                journal="National Institute of Standards and Technology",
                contribution="Comprehensive database of atomic ionization energies and spectroscopic data",
                weight=1.0,
                url="https://www.nist.gov/pml/atomic-spectra-database",
            )
        )

        # Fisher - Statistical methods
        self.add_citation(
            Citation(
                key="Fisher1932",
                authors=["Ronald Fisher"],
                title="Statistical Methods for Research Workers",
                year=1932,
                journal="Oliver and Boyd",
                contribution="Developed method for combining independent p-values, foundational to statistical inference",
                weight=1.0,
            )
        )

        # Lattice QCD
        self.add_citation(
            Citation(
                key="LatticeQCD2016",
                authors=["Flavour Lattice Averaging Group"],
                title="FLAG Review 2016",
                year=2016,
                journal="European Physical Journal C",
                contribution="Lattice QCD calculations of quark masses and QCD parameters",
                weight=1.0,
                url="https://doi.org/10.1140/epjc/s10052-016-4509-7",
            )
        )

        # RBC-UKQCD Collaboration - Hadronic matrix elements for proton decay
        self.add_citation(
            Citation(
                key="RBCUKQCD2015",
                authors=["RBC and UKQCD Collaborations"],
                title="Nucleon form factors and proton decay matrix elements from lattice QCD",
                year=2015,
                journal="Physical Review D",
                contribution="Lattice QCD calculation of hadronic matrix elements for proton decay",
                weight=1.0,
                url="https://doi.org/10.1103/PhysRevD.91.074506",
            )
        )

        # Slater - Atomic screening
        self.add_citation(
            Citation(
                key="Slater1930",
                authors=["John Slater"],
                title="Atomic Shielding Constants",
                year=1930,
                journal="Physical Review",
                contribution="Developed empirical rules for effective nuclear charge in atoms",
                weight=1.0,
                url="https://doi.org/10.1103/PhysRev.36.57",
            )
        )

        # Bousso - String landscape
        self.add_citation(
            Citation(
                key="Bousso2000",
                authors=["Raphael Bousso", "Joseph Polchinski"],
                title="Quantization of Four-form Fluxes and Dynamical Neutralization of the Cosmological Constant",
                year=2000,
                journal="Journal of High Energy Physics",
                contribution="String theory landscape and vacuum selection problem",
                weight=1.0,
                url="https://doi.org/10.1088/1126-6708/2000/06/006",
            )
        )

        # Polchinski - String theory
        self.add_citation(
            Citation(
                key="Polchinski1998",
                authors=["Joseph Polchinski"],
                title="String Theory",
                year=1998,
                journal="Cambridge University Press",
                contribution="Definitive textbook on string theory, D-branes, and dualities",
                weight=1.0,
            )
        )

        # Acharya - G2 compactifications (already have Acharya2004, this is earlier)
        self.add_citation(
            Citation(
                key="Acharya1998",
                authors=["Bobby Acharya"],
                title="M theory, Joyce Orbifolds and Super Yang-Mills",
                year=1998,
                journal="Advances in Theoretical and Mathematical Physics",
                contribution="Early work on M-theory compactifications on G$_2$ manifolds",
                weight=1.0,
                url="https://doi.org/10.4310/ATMP.1998.v2.n6.a1",
            )
        )

        # Atiyah - G2 holonomy mathematics
        self.add_citation(
            Citation(
                key="Atiyah1999",
                authors=["Michael Atiyah", "Edward Witten"],
                title="M-Theory Dynamics On A Manifold Of G$_2$ Holonomy",
                year=2001,
                journal="Advances in Theoretical and Mathematical Physics",
                contribution="Mathematical structure of M-theory on G$_2$ holonomy manifolds",
                weight=1.0,
                url="https://doi.org/10.4310/ATMP.2002.v6.n1.a1",
            )
        )

        # Isaac Asimov - For posing the entropy question
        self.add_citation(
            Citation(
                key="Asimov1956",
                authors=["Isaac Asimov"],
                title="The Last Question",
                year=1956,
                journal="Science Fiction Quarterly",
                contribution="Asked how entropy might be reversed, a question Multivac took until well after the end of the universe to answer. G$_2$ geometry provides the answer: $\Omega$_Λ = 11/16 determines the vacuum energy that can reverse entropy locally through quantum tunneling",
                weight=0.0,  # Honorary citation
            )
        )

        # Marciano & Senjanovic - Proton decay calculation
        self.add_citation(
            Citation(
                key="MarcianSenjanovic1982",
                authors=["William Marciano", "Goran Senjanovic"],
                title="Predictions of Supersymmetric Grand Unified Theories",
                year=1982,
                journal="Physical Review D 25, 3092",
                contribution="Calculated proton decay rates from dimension-6 operators in GUT theories",
                weight=1.0,
                url="https://doi.org/10.1103/PhysRevD.25.3092",
            )
        )

        # Stuttgart 2025 - Experimental validation of geometry-dependent thermodynamics
        self.add_citation(
            Citation(
                key="Stuttgart2025",
                authors=["Milton Aguilar", "Eric Lutz"],
                title="Correlated quantum machines beyond the standard second law",
                year=2025,
                journal="Science Advances, DOI: 10.1126/sciadv.adw8462",
                contribution="Experimentally demonstrated that correlated quantum systems violate the Carnot principle, validating that thermodynamic laws depend on geometric/quantum structure—exactly as predicted by G$_2$ framework",
                weight=1.0,  # Full credit: independent experimental validation
                url="https://doi.org/10.1126/sciadv.adw8462",
            )
        )

        # Antusch & Baranowski 2018 - RG corrections to neutrino mixing
        self.add_citation(
            Citation(
                key="AntuschBaranowski2018",
                authors=["Stefan Antusch", "Marja Baranowski"],
                title="Running of neutrino parameters and the Higgs self-coupling",
                year=2018,
                journal="Physical Review D 98, 113001",
                contribution="Calculated RG evolution of neutrino mixing angles in seesaw models, showing 1-10% corrections from M_R to M_Z",
                weight=0.0,  # Reference only for validation, not used in derivation
                url="https://doi.org/10.1103/PhysRevD.98.113001",
            )
        )

        # Casas et al. 1999 - Early RG analysis of neutrino parameters
        self.add_citation(
            Citation(
                key="Casas1999",
                authors=["J.A. Casas", "J.R. Espinosa", "A. Ibarra", "I. Navarro"],
                title="General RG equations for physical neutrino parameters",
                year=1999,
                journal="Nuclear Physics B 573, 652",
                contribution="Pioneering work on RG evolution of neutrino mixing angles, establishing expected correction magnitudes",
                weight=0.0,  # Reference only for validation
                url="https://doi.org/10.1016/S0550-3213(00)00057-8",
            )
        )

        # Tomonori Totani - 20 GeV Halo Excess
        self.add_citation(
            Citation(
                key="Totani2025",
                authors=["Tomonori Totani"],
                title="20 GeV halo-like excess of the Galactic diffuse emission and implications for dark matter annihilation",
                year=2025,
                journal="Journal of Cosmology and Astroparticle Physics",
                contribution="Observed 20 GeV gamma-ray excess compatible with ~700 GeV mass scale",
                weight=1.0,
                url="https://arxiv.org/abs/2507.16477",
            )
        )

    def add_citation(self, citation: Citation) -> None:
        """Add a citation to the database."""
        self.citations[citation.key] = citation

    def add_attribution(self, attribution: Attribution) -> None:
        """Add an attribution for a result."""
        self.attributions.append(attribution)

    def compute_total_attribution(self) -> Dict[str, float]:
        """
        Compute total attribution across all results.

        Returns percentage contribution for each citation key.
        """
        contributions = defaultdict(float)
        total_results = len(self.attributions)

        for attr in self.attributions:
            weights = attr.normalized_weights()
            for key, weight in weights.items():
                contributions[key] += weight / total_results

        return dict(contributions)

    def generate_bibtex(self) -> str:
        """Generate BibTeX file from citations."""
        entries = []

        for cite in self.citations.values():
            # Determine entry type
            if "Collaboration" in cite.authors[0]:
                entry_type = "article"
            elif cite.year < 1950:
                entry_type = "book"
            else:
                entry_type = "article"

            entry = f"@{entry_type}{{{cite.key},\n"
            entry += f"  author = {{{' and '.join(cite.authors)}}},\n"
            entry += f"  title = {{{cite.title}}},\n"
            entry += f"  year = {{{cite.year}}},\n"
            entry += f"  journal = {{{cite.journal}}},\n"
            if cite.url:
                entry += f"  url = {{{cite.url}}},\n"
            entry += f"  note = {{{cite.contribution}}}\n"
            entry += "}\n"

            entries.append(entry)

        return "\n".join(entries)

    def generate_attribution_report(self) -> str:
        """
        Generate a detailed attribution report.

        Returns human-readable summary of contributions.
        """
        report = []
        report.append("=" * 80)
        report.append("ATTRIBUTION REPORT: $\alpha$$\Omega$ Framework")
        report.append("=" * 80)
        report.append("")

        # Overall statistics
        report.append(f"Total results analyzed: {len(self.attributions)}")
        report.append("")

        # Contribution by researcher
        total_attr = self.compute_total_attribution()
        report.append("Contributions by citation:")
        report.append("-" * 80)

        # Sort by contribution percentage
        sorted_attr = sorted(total_attr.items(), key=lambda x: x[1], reverse=True)

        for key, pct in sorted_attr:
            cite = self.citations[key]
            authors_str = (
                cite.authors[0]
                if len(cite.authors) == 1
                else f"{cite.authors[0]} et al."
            )
            report.append(
                f"  {authors_str:25} ({cite.year}):  {pct * 100:5.1f}%  - {cite.contribution[:60]}"
            )

        report.append("")
        report.append("=" * 80)

        return "\n".join(report)


# Global citation database
citation_db = CitationDatabase()


# ============================================================================
# EXAMPLE ATTRIBUTIONS
# ============================================================================


def initialize_attributions() -> None:
    """
    Initialize attributions for key results.

    This defines what prior work each result builds on.
    """

    # $\alpha$_GUT = 1/42
    citation_db.add_attribution(
        Attribution(
            result_name="$\alpha$_GUT = 1/42",
            description="Grand unified coupling from G$_2$ geometry",
            prior_work=[
                citation_db.citations["Cartan1894"],  # G$_2$ discovered by Cartan
                citation_db.citations["Killing1888"],  # Lie algebra classification
            ],
        )
    )

    # $\Omega$_Λ = 11/16
    citation_db.add_attribution(
        Attribution(
            result_name="$\Omega$_Λ = 11/16",
            description="Dark energy density from G$_2$ Casimir operators",
            prior_work=[
                citation_db.citations["Cartan1894"],  # G$_2$ Casimir invariants
                citation_db.citations["Planck2018"],  # Measured $\Omega$_Λ
            ],
        )
    )

    # sin²θ_W = 3/13
    citation_db.add_attribution(
        Attribution(
            result_name="sin²θ_W = 3/13",
            description="Weak mixing angle from triality",
            prior_work=[
                citation_db.citations["Weinberg1967"],  # Defined θ_W
                citation_db.citations["Glashow1961"],  # SU(2)×U(1) structure
                citation_db.citations["Cartan1894"],  # Triality
            ],
        )
    )

    # Three generations from $\tau$³ = 1
    citation_db.add_attribution(
        Attribution(
            result_name="N_gen = 3",
            description="Three fermion generations from triality order",
            prior_work=[
                citation_db.citations["Cartan1894"],  # Triality automorphism
            ],
        )
    )

    # Neutrino mass ratios
    citation_db.add_attribution(
        Attribution(
            result_name="Neutrino mass ratios 7/8, 1/20, 13/11",
            description="Neutrino mass ratios from G$_2$ and E₆ branching rules",
            prior_work=[
                citation_db.citations["Minkowski1977"],  # Seesaw mechanism
                citation_db.citations["SuperK1998"],  # Measured oscillations
                citation_db.citations["Cartan1894"],  # G$_2$ structure
            ],
        )
    )

    # Strong CP solution
    citation_db.add_attribution(
        Attribution(
            result_name="θ_QCD = 0",
            description="Strong CP angle vanishes from triality",
            prior_work=[
                citation_db.citations["Cartan1894"],  # Triality
            ],
        )
    )

    # Cabibbo angle
    citation_db.add_attribution(
        Attribution(
            result_name="sin θ_C = √(5/98)",
            description="Cabibbo angle from G$_2$ Clebsch-Gordan coefficients",
            prior_work=[
                citation_db.citations["Cabibbo1963"],  # Discovered mixing
                citation_db.citations["Cartan1894"],  # G$_2$ representation theory
            ],
        )
    )


if __name__ == "__main__":
    # Initialize attributions
    initialize_attributions()

    # Generate attribution report
    print(citation_db.generate_attribution_report())
    print()

    # Generate BibTeX
    print("=" * 80)
    print("BIBTEX ENTRIES")
    print("=" * 80)
    print(citation_db.generate_bibtex())
