# Planning without an ordinary power calculation

Some methods are sized by precision, some study sizes are fixed, some rules only provide a model
capacity floor, and some designs require a bespoke simulation. Other procedures, such as an
imputation or sensitivity analysis, are not sized independently at all.

The two scripts make those alternatives executable and explicit. A caller supplies the analysis
entry's reason and the route catalogue's required inputs and next step. They deliberately return no
sample size. A plan that needs bespoke simulation must specify the data-generating model and the
analysis fitted inside every replicate, then report the seed, attempted replicates, successful fits,
failed fits, estimated power and Monte Carlo standard error.
