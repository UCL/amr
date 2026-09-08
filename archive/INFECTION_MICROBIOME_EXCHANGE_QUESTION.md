# Clarification needed: infection–microbiome exchange

Checked against **official UCL/amr `main`**, commit `a0f62f699cd2e178c7359659c49b4dadfc5c2a42` (7 September 2026), on 8 September 2026.

## Confirmed behaviour

The model contains an operation for exchanging resistance characteristics between an infection and carriage of the same organism.

In the checked code, the exchange block requires a **positive infection level**, but sits inside the branch for **no positive infection episode**. Under that control flow, the exchange does not execute.

The current model description, section 8.2, already acknowledges that this pathway is inactive and that changing its transfer-probability parameter cannot activate it. This is therefore **not a newly discovered, unacknowledged defect report**.

Earlier documentation described daily bidirectional exchange when the same organism was present in both compartments. The latest documentation describes the inactive implementation without changing the exchange block. That establishes current behaviour, but does not settle whether inactivity is intentional or an implementation limitation.

## Response requested

1. **Should this exchange remain inactive, or should it operate when infection and carriage coexist?**
2. If it should operate, is it an explicit exception to the rule that a positive infection episode freezes carriage acquisition, clearance, resistance progression and scheduling?

Please identify the intended rule and any existing decision or regression test that establishes it. Carrier-to-new-infection inheritance and inter-species horizontal transfer are separate pathways; this question is only about the same-organism infection–carriage exchange.

No parameter difference is being reported as a defect. The aggregate effect of enabling this pathway has not been measured by this inspection.

## Original-repository source locations

- [Outer episode-ownership branch](https://github.com/UCL/amr/blob/a0f62f699cd2e178c7359659c49b4dadfc5c2a42/src/rules/mod.rs#L5201)
- [Exchange block](https://github.com/UCL/amr/blob/a0f62f699cd2e178c7359659c49b4dadfc5c2a42/src/rules/mod.rs#L5526)
- [Current model description, section 8.2](https://github.com/UCL/amr/blob/a0f62f699cd2e178c7359659c49b4dadfc5c2a42/model_description/MODEL_DESCRIPTION.md#L1875)
