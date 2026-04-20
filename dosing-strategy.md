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

# Aminoglycoside Extended-Interval (EI) Dosing Algorithm

```mermaid
flowchart TD
    %% Step 1: Checklist
    Start[<b>Step 1: EI Eligibility Checklist</b><br/>- Gram-positive synergy<br/>- Renal insufficiency / AKI / CrCl 20<br/>- Hemodialysis / CRRT<br/>- Surgical prophylaxis<br/>- Pregnancy / Neonatal<br/>- NTM infection] --> B{Any Exclusions<br/>Present?}
    
    %% Step 2: Outcomes
    B -- "No" --> EI[<b>Extended-Interval Dosing</b>]
    B -- "Yes" --> C{Is it Synergy, AKI,<br/>or CrCl 20?}
    
    C -- "Yes" --> Conv[<b>Conventional Dosing</b>]

    %% Step 3: Referral
    C -- "No" --> D[<b>Refer to Specific Dosing Section</b><br/>or Contact ID Pharmacy<br/><br/>- Hemodialysis<br/>- CRRT<br/>- Surgical Prophylaxis<br/>- Neonatal Population<br/>- NTM Infections]

    %% STYLING SECTION
    %% This 'class' approach is the most stable way to force width
    classDef wideBox width:400px,text-align:left;
    classDef outcomeBox width:300px,text-align:center;

    class Start,D wideBox;
    class EI,Conv outcomeBox;

    %% Individual Colors
    style Start fill:#fff9c4,stroke:#fbc02d
    style EI fill:#d4edda,stroke:#28a745,stroke-width:2px
    style Conv fill:#f8d7da,stroke:#dc3545,stroke-width:2px
    style D fill:#e1f5fe,stroke:#01579b
```
