# Aminoglycoside Extended-Interval (EI) Dosing Algorithm

---

## 1. Eligibility for Extended-Interval (EI) Dosing

Patients are eligible for EI dosing unless any exclusion criteria below are present.

---

## 2. Exclusion Criteria for EI Dosing

- Gram-positive synergy indication (e.g., endocarditis)
- Renal insufficiency requiring conventional dosing consideration
- Hemodialysis (HD)
- Continuous renal replacement therapy (CRRT)
- Acute kidney injury with unstable renal function
- Surgical prophylaxis
- Pregnancy or neonatal population
- Nontuberculous mycobacterial (NTM) infection

---

## 3. Dosing Pathways

- **Gram-positive synergy (e.g., endocarditis)** → Conventional dosing
- **Renal insufficiency / AKI** → Conventional dosing
- **Hemodialysis (HD)** → Refer to HD dosing protocol
- **CRRT** → Refer to CRRT dosing protocol
- **Surgical prophylaxis** → Refer to surgical prophylaxis guideline
- **Neonatal population** → Refer to neonatal dosing section
- **NTM infections** → Refer to NTM dosing section

---

## 4. Flowchart
# Aminoglycoside Dosing Algorithm

```mermaid
flowchart TD
    %% Step 1: Initial Triage
    Triage{"Step 1: Primary Triage
    Adult Patient AND
    Stable Renal Function AND
    CrCl > 20 mL/min?"} 
    
    %% Standard Path
    Triage -- "Yes" --> Start[Step 2: Check Final Exclusions
    - Gram-positive synergy
    - Renal Replacement Therapy
    - HD / PD / CRRT
    - Surgical Prophylaxis
    - Neonatal Population
    - NTM / Mycobacterial Infections]

    %% Immediate Diversion
    Triage -- "No" --> C

    %% Final Check for EI
    Start --> B{Any Exclusions
    Present?}
    
    B -- No --> EI[Extended-Interval Dosing]
    B -- Yes --> C{Is it Synergy, AKI,
    or CrCl < 20?}
    
    %% Outcomes
    C -- Yes --> Conv[Conventional Dosing]
    
    C -- No --> D[Refer to Specific Dosing Section
    or Contact ID Pharmacy
    
    - Hemodialysis / PD
    - CRRT
    - Surgical Prophylaxis
    - Neonatal Population
    - NTM Infections]

    %% STYLING SECTION
    classDef wideBox min-width:400px,text-align:left;
    classDef outcomeBox min-width:300px,text-align:center;
    classDef diamond padding:20px;

    class Start,D wideBox;
    class EI,Conv outcomeBox;
    class Triage,B,C diamond;

    %% Individual Colors
    style Triage fill:#e8eaf6,stroke:#3f51b5
    style Start fill:#fff9c4,stroke:#fbc02d
    style EI fill:#d4edda,stroke:#28a745,stroke-width:2px
    style Conv fill:#f8d7da,stroke:#dc3545,stroke-width:2px
    style D fill:#e1f5fe,stroke:#01579b
```
