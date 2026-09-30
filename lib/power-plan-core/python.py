"""Standard power and precision planning selected by a typed calc key.

Fixture values are synthetic. Generated code must replace them with sourced study assumptions.
Rule-of-thumb floors remain explicitly distinct from power calculations.
"""

import csv
import math
from scipy.stats import f, ncf, nct, norm, t


def number(value, default=math.nan):
    try:
        if value is None or value == "":
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def two_sided_t_power(ncp_value, df, alpha):
    critical = t.ppf(1 - alpha / 2, df)
    return nct.cdf(-critical, df, ncp_value) + nct.sf(critical, df, ncp_value)


def exact_t_n(delta, sd, alpha, target, ratio=1, two_groups=False):
    z = norm.ppf(1 - alpha / 2) + norm.ppf(target)
    if two_groups:
        start = max(2, math.floor((1 + 1 / ratio) * z**2 * sd**2 / delta**2) - 3)
        for n1 in range(start, 100001):
            n2 = max(2, math.ceil(ratio * n1))
            se = sd * math.sqrt(1 / n1 + 1 / n2)
            if two_sided_t_power(abs(delta) / se, n1 + n2 - 2, alpha) >= target:
                return {"n1": n1, "n2": n2, "total": n1 + n2}
    else:
        start = max(2, math.floor(z**2 * sd**2 / delta**2) - 3)
        for n in range(start, 100001):
            if two_sided_t_power(abs(delta) * math.sqrt(n) / sd, n - 1, alpha) >= target:
                return {"total": n}
    raise ValueError("target power was not reached within the search bound")


def power_plan(calc, x):
    alpha = number(x.get("alpha"), 0.05)
    target = number(x.get("power"), 0.8)
    ratio = number(x.get("ratio"), 1)
    za, zb = norm.ppf(1 - alpha / 2), norm.ppf(target)
    if calc == "twoMeans":
        return exact_t_n(number(x["delta"]), number(x["sd"]), alpha, target, ratio, True)
    if calc in {"paired", "oneMean"}:
        return exact_t_n(number(x["delta"]), number(x["sd"]), alpha, target)
    if calc in {"twoProportions", "oneProportion"}:
        p1, p2 = number(x["p1"]), number(x["p2"])
        difference = abs(p1 - p2)
        if calc == "twoProportions":
            pooled = (p1 + ratio * p2) / (1 + ratio)
            a = za * math.sqrt((1 + 1 / ratio) * pooled * (1 - pooled))
            b = zb * math.sqrt(p1 * (1 - p1) + p2 * (1 - p2) / ratio)
            exact = (a + b) ** 2 / difference**2
            n1, n2 = math.ceil(exact), math.ceil(ratio * exact)
            return {"n1": n1, "n2": n2, "total": n1 + n2}
        exact = (za * math.sqrt(p1 * (1 - p1)) + zb * math.sqrt(p2 * (1 - p2))) ** 2 / difference**2
        return {"total": math.ceil(exact)}
    if calc == "correlation":
        fisher = math.atanh(abs(number(x["r"])))
        return {"total": math.ceil(((za + zb) / fisher) ** 2 + 3)}
    if calc == "logrank":
        hazard_ratio, event_rate = number(x["hr"]), number(x["eventRate"])
        pi1 = 1 / (1 + ratio)
        events = (za + zb) ** 2 / (pi1 * (1 - pi1) * math.log(hazard_ratio) ** 2)
        return {"events": math.ceil(events), "total": math.ceil(events / event_rate)}
    if calc == "logisticEPV":
        events = 10 * round(number(x["predictors"]))
        return {"events": events, "total": math.ceil(events / number(x["eventRate"]))}
    if calc == "linearSPV":
        predictors = round(number(x["predictors"]))
        return {"total": max(50 + 8 * predictors, 104 + predictors)}
    if calc == "dxAccuracy":
        sensitivity, specificity = number(x["sens"]), number(x["spec"])
        prevalence, width = number(x["prev"]), number(x["width"])
        diseased = math.ceil(za**2 * sensitivity * (1 - sensitivity) / width**2)
        non_diseased = math.ceil(za**2 * specificity * (1 - specificity) / width**2)
        return {"diseased": diseased, "nonDiseased": non_diseased,
                "total": max(math.ceil(diseased / prevalence), math.ceil(non_diseased / (1 - prevalence)))}
    if calc == "agreementKappa":
        kappa, prevalence, width = number(x["kappa"]), number(x["prev"]), number(x["width"])
        variance = (1 - kappa) * ((1 - kappa) * (1 - 2 * kappa) + kappa * (2 - kappa) / (2 * prevalence * (1 - prevalence)))
        return {"total": max(20, math.ceil(za**2 * variance / width**2))}
    if calc == "anovaKGroups":
        groups, effect = round(number(x["groups"])), number(x["f"])
        for per_group in range(2, 100001):
            total = groups * per_group
            critical = f.ppf(1 - alpha, groups - 1, total - groups)
            attained = ncf.sf(critical, groups - 1, total - groups, effect**2 * total)
            if attained >= target:
                return {"perGroup": per_group, "total": total}
    if calc == "proportionPrecision":
        p, width, icc = number(x["p"]), number(x["width"]), number(x.get("icc"), 0)
        clusters, cluster_size = number(x.get("clusters")), number(x.get("clusterSize"))
        effective = za**2 * p * (1 - p) / width**2
        if math.isfinite(cluster_size):
            exact = effective * (1 + (cluster_size - 1) * icc)
        elif math.isfinite(clusters):
            denominator = 1 - effective * icc / clusters
            if denominator <= 0:
                raise ValueError("more clusters are required to reach this precision")
            exact = effective * (1 - icc) / denominator
        else:
            exact = effective
        total = math.ceil(exact - 1e-9)
        if math.isfinite(clusters):
            total = max(total, round(clusters))
        return {"total": total}
    if calc == "noninferiorityProportions":
        p_ref = number(x["pRef"]); p_new = number(x.get("pNew"), p_ref)
        margin = number(x["margin"]); one_alpha = number(x.get("alpha"), 0.025)
        expected = p_ref - p_new if str(x.get("worseIs", "higher")).lower() == "lower" else p_new - p_ref
        room = margin - expected
        n_ref = math.ceil((norm.ppf(1 - one_alpha) + norm.ppf(target)) ** 2 *
                          (p_ref * (1 - p_ref) + p_new * (1 - p_new) / ratio) / room**2)
        n_new = math.ceil(n_ref * ratio)
        return {"nReference": n_ref, "nNew": n_new, "total": n_ref + n_new}
    raise ValueError(f"unsupported calc {calc}")


with open("fixture.csv", newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle))
inputs_by_calc = {}
for row in rows:
    inputs_by_calc.setdefault(row["calc"], {})[row["key"]] = row["value"]
results = {}
for calc, inputs in inputs_by_calc.items():
    for name, value in power_plan(calc, inputs).items():
        results[f"{calc}_{name}"] = value

print("These are analyzable sample sizes before attrition or missing-data inflation.")
print("logisticEPV and linearSPV are rule-of-thumb floors, not power calculations.")
print("\n--- HARNESS ---")
for name, value in results.items():
    print(f"{name}={value:.10f}")
