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

```mermaid
flowchart TD
    %% Step 1: Initial Triage
    Triage{"<b>Step 1: Primary Triage</b><br/>Adult Patient AND<br/>Stable Renal Function AND<br/>CrCl > 20 mL/min?"} 
    
    %% Standard Path
    Triage -- "Yes" --> Start[<b>Step 2: Check Final Exclusions</b><br/>- Gram-positive synergy<br/>- Renal Replacement Therapy / HD / CRRT<br/>- Surgical Prophylaxis<br/>- Neonatal Population<br/>- NTM / Mycobacterial Infections]

    %% Immediate Diversion for non-standard adults
    Triage -- "No" --> C

    %% Final Check for EI
    Start --> B{Any Exclusions<br/>Present?}
    
    B -- No --> EI[<b>Extended-Interval Dosing</b>]
    B -- Yes --> C{Is it Synergy, AKI,<br/>or CrCl < 20?}
    
    %% Outcomes
    C -- Yes --> Conv[<b>Conventional Dosing</b>]
    
    C -- No --> D[<b>Refer to Specific Dosing Section</b><br/>or Contact ID Pharmacy<br/><br/>- Renal Replacement (HD/PD/CRRT)<br/>- Surgical Prophylaxis<br/>- Neonatal Population<br/>- NTM Infections]

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
