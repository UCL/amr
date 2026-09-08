I have implemented the agreed infection–carriage exchange change. At the start of each organism's daily update, an existing positive infection episode, including a fading episode, and separate carriage of that organism receive at most one exchange draw when their eligible profiles differ. It uses the existing daily probability and counterfactual scaling; a zero multiplier disables it.

A successful exchange shares their existing eligible characteristics and raises the corresponding resistance measures where needed. It preserves compartment ownership and episode dates and does not itself promote the predominant-strain profile. New infections use the separate inheritance route on their acquisition day and first become eligible for exchange the following day.

Ordinary carriage acquisition, clearance, emergence and reversion remain suspended during infection. Changing that broader rule remains a separate decision. The model description and checklist now reflect this scope.

All eight exchange tests pass. Full verification reports 215 Rust and 207 Python passes, with the same two existing Rust parameter-assertion failures and one Python whitespace-assertion failure.

This changes model dynamics without retuning the exchange parameter. Its calibration impact remains unmeasured and requires new simulation runs.
