# Stochastic Simulation of Histone-Based Epigenetic Memory
### A Study Inspired by Dodd et al. (Cell, 2007)

---

## 📌 Project Overview

This project presents a **stochastic Monte Carlo simulation** of histone-based epigenetic memory using a theoretical model proposed by **Dodd et al. (Cell, 2007)**. The model explores how nucleosome modification states — Methylated (M), Unmodified (U), and Acetylated (A) — can give rise to **bistability** and **stable epigenetic memory** through feedback-driven recruitment and random noise.

The simulation was independently implemented in Python as part of an **On-Job Training (OJT)** at the **Physical Biology Lab, IIT Bombay**, under the guidance of **Dr. Ranjith Padinhateeri**.

---

## 🏛️ Institutional Details

| | |
|---|---|
| **Institution** | SIES College of Arts, Science & Commerce (Autonomous), Sion (W), Mumbai |
| **Training Lab** | Physical Biology Lab, IIT Bombay |
| **Supervisor** | Dr. Ranjith Padinhateeri, Department of Biosciences & Bioengineering, IIT Bombay |
| **Mentor** | Mr. Vinoth M., Ph.D. Scholar, Physical Biology Lab |
| **Authors** | Yadav Janhavi Tarendrakumar & Chowdhry Nida Abu Bilal |
| **Program** | M.Sc. Physics |
| **Year** | 2025 |

---

## 🧬 Biological Background

Epigenetic memory refers to the ability of cells to remember and maintain specific patterns of gene expression across cell divisions — **without any change in the DNA sequence**.

A key mechanism behind this involves **nucleosome modifications**:
- **Methylation (M)** → gene silencing (e.g., H3K9me)
- **Acetylation (A)** → gene activation (e.g., H3K14Ac)
- **Unmodified (U)** → neutral/intermediate state

Modified nucleosomes recruit enzymes that reinforce the **same modification in neighboring nucleosomes**, forming **positive feedback loops** — the core of epigenetic bistability.

---

## 🎯 Objectives

1. Recreate the stochastic nucleosome modification model from Dodd et al. (2007) in Python
2. Study the effect of **feedback-to-noise ratio (F)** on bistability
3. Measure **gap score (G)** and **state lifetime** as metrics of epigenetic memory
4. Analyze the role of **cooperativity** in bistability
5. Study the effect of **spatial constraints** on memory formation
6. Compare simulation results with the original paper

---

## 🧪 The Model

- **System**: Linear array of N = 60 nucleosomes
- **States**: Each nucleosome ∈ {M, U, A}
- **Transitions**: Only via U (no direct M ↔ A)
- **Feedback**: Modified nucleosomes recruit enzymes to convert neighbors
- **Noise**: Random spontaneous state changes
- **Key Parameter**: Feedback-to-noise ratio F = α / (1 − α)

---

## 📂 Project Structure

```
epigenetic-memory-dodd2007/
│
├── simulation/
│   ├── bistability.py                  # Time evolution & bistability (Section 5.1)
│   ├── lifetime_gap_score.py           # Lifetime & Gap Score vs F (Section 5.2 & 5.3)
│   ├── cooperativity_gap_score.py      # G vs F for Cases A, B, C (Section 5.4)
│   ├── cooperativity_distribution.py   # P(M-A) distribution under cooperativity (Section 5.4)
│   ├── spatial_constraints.py          # G vs F for spatial models (Section 5.5)
│   └── spatial_distribution.py         # P(M-A) under spatial constraints (Section 5.5)
│
├── results/
│   └── (simulation graphs and figures)
│
├── requirements.txt
└── README.md
```

---

## 📊 Simulations Performed

### 1. Bistability vs Feedback-to-Noise Ratio
- Simulated nucleosome dynamics at F = 0.4, 1.0, 1.4, 2.0
- Plotted time traces of M(t) and probability distributions P(M)
- **Key finding**: Clear bimodal distribution emerges at high F → strong bistability

### 2. Lifetime vs Feedback-to-Noise Ratio
- Measured average duration of high-M and high-A states
- Used 1.5× threshold rule to classify dominant states
- **Key finding**: Lifetime increases exponentially with F

### 3. Gap Score vs Feedback-to-Noise Ratio
- Gap Score G = |M − A| / (M + A)
- Sigmoidal G vs F curve with sharp transition around F = 1.0–1.5
- **Key finding**: G → 1 at high F, confirming robust bistability

### 4. Effect of Cooperativity (Cases A, B, C)
| Case | Feedback Type | Result |
|------|--------------|--------|
| A | Both modification + demodification | Strong bistability even without explicit cooperativity |
| B | Modification only | Cooperativity essential for bistability |
| C | Demodification only | Even cooperativity insufficient for strong bistability |

### 5. Spatial Constraint Effects (Cases A, B, C)
| Case | Spatial Model | Result |
|------|--------------|--------|
| A | No constraint (global) | Strong bistability |
| B | Neighbor-only recruitment | Weak bistability, slow G rise |
| C | Power-law decay (∝ 1/d^1.5) | Moderate bistability restored |

---

## 🛠️ Libraries Used

| Library | Purpose |
|---------|---------|
| `numpy` | Array operations, random number generation |
| `matplotlib` | Plotting graphs and distributions |
| `collections` | Counter for histogram building |

---

## ▶️ How to Run

**Step 1: Install dependencies**
```bash
pip install numpy matplotlib
```

**Step 2: Run any simulation**
```bash
python simulation/bistability.py
python simulation/lifetime_gap_score.py
python simulation/cooperativity_gap_score.py
python simulation/cooperativity_distribution.py
python simulation/spatial_constraints.py
python simulation/spatial_distribution.py
```

> ⚠️ Note: Some simulations (especially cooperativity and spatial) are computationally heavy and may take several minutes to hours depending on your system.

---

## 📈 Key Results

- Bistability **requires** a minimum feedback-to-noise ratio (F ≥ ~1.5)
- **Full feedback model** (Case A) shows bistability even without explicit cooperativity due to implicit two-step recruitment
- **Modification-only feedback** (Case B) requires explicit cooperativity for stable memory
- **Neighbor-only spatial constraints** significantly reduce bistability — long-range interactions are crucial
- **Power-law recruitment** partially restores memory, consistent with 3D chromatin folding

---

## 📖 Reference

> Dodd, I. B., Micheelsen, M. A., Sneppen, K., & Thon, G. (2007).
> **Theoretical Analysis of Epigenetic Cell Memory by Nucleosome Modification.**
> *Cell*, 129(4), 813–822.
> https://doi.org/10.1016/j.cell.2007.02.053

---

## 🙏 Acknowledgements

We thank **Dr. Ranjith Padinhateeri** (IIT Bombay) for providing the opportunity to work in his lab and for his inspiring guidance. We are grateful to **Mr. Vinoth M.** (Ph.D. Scholar) for his constant mentorship and support throughout this project.

We also thank **SIES College of Arts, Science & Commerce (Autonomous)** for facilitating this OJT opportunity.

---

*M.Sc. Physics | SIES College, Mumbai | OJT at Physical Biology Lab, IIT Bombay | 2025*
