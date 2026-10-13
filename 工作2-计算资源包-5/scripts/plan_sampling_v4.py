"""Reproduce the advisory, structure-informed v4 microsolvation sampling plan.

This writes planning metadata only.  It neither launches GOAT nor asserts that
any proposed number of starts or minima establishes conformational convergence.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


SCENARIO = "DMSO_02_W01_trace_water_sensitivity"
PLAN_VERSION = "2026-10-14-v4"

# Counts are deliberate per-structure decisions, not a rotatable-bond formula.
# Ring conformational freedom and solvent-binding alternatives matter even when
# RDKit reports no acyclic rotatable bonds.
DECISIONS = {
    "F2_E01": (8, "苯乙烯环氧化物较刚性；采样重点为环氧氧周围不同溶剂取向。"),
    "F2_E03": (8, "对氯取代不增加主骨架自由度；保留环氧端与芳环面的溶剂起点。"),
    "F2_E05": (12, "腈基和环氧氧提供不同水结合位置，较单一环氧化物增加起点。"),
    "F2_E09": (16, "正己基链有5个可旋转键；需要链构象与环氧端溶剂布置的交叉采样。"),
    "F2_E11": (8, "桥环骨架较刚性；endo/exo分别各8个起点，不互相混合。"),
    "F2_E13": (16, "八元环有明显环构象自由度；不能因RDKit可旋转键为0而少采样。"),
    "F2_E12": (12, "羟基和环氧氧存在竞争氢键及水桥；增加局部溶剂布局。"),
    "F2_E14": (16, "丙基侧链、羟甲基和双氧受体形成多个构象/氢键组合。"),
    "F2_E15": (16, "苯基取向、羟甲基与环氧氧可竞争溶剂结合。"),
    "F2_E17": (24, "柔性萜烯侧链与羟甲基并存；需要折叠/伸展和不同氢键构型。"),
    "F2_E19": (16, "环己烷椅式及羟基/环氧端取向均需覆盖。"),
    "F2_E21": (16, "稠合环构象加羟乙基侧链；重点比较溶剂桥与独立结合。"),
    "F2_E22": (12, "环酮、环氧氧和环构象均需采样；两个面异构模型各12个独立起点。"),
    "F2_E23": (16, "两个苯基取向、羰基和环氧端组合较多。"),
    "F2_E24": (12, "稠环较刚性，但两个羰基与环氧氧提供三个主要受体区域。"),
    "F2_E18": (32, "9个可旋转键、长多烯链、酯和醇；本组最高柔性层级。"),
    "F2_E20": (32, "8个可旋转键及两个羟基；链折叠和内/外氢键竞争显著。"),
    "F2_S01": (8, "单芳基亚砜主骨架简单；重点覆盖S=O周围不同氢键取向。"),
    "F2_S06": (12, "两个苯基绕硫旋转和S=O暴露面的组合。"),
    "F2_S07": (12, "苄基引入额外扭转；溶剂结合与苯环取向共同采样。"),
    "F2_S19": (12, "硝基与亚砜氧构成不同溶剂结合区域。"),
    "F2_S22": (12, "环丙基/芳基相对取向及S=O配位面需要分别覆盖。"),
    "F2_S23": (12, "乙烯基/芳基亚砜扭转与溶剂取向需要交叉采样。"),
    "F2_S12": (24, "两条丁基链共6个可旋转键，且砜有两个受体氧。"),
    "F2_S24": (16, "烯丙基链、芳基扭转及两个砜氧结合位置并存。"),
    "F2_S25": (16, "羟乙基与砜氧可形成内氢键或水桥；需覆盖竞争模式。"),
    "F2_S11": (8, "稠合芳香骨架刚性，主要变化为两个砜氧的溶剂位置。"),
    "F2_S26": (12, "噻蒽中心环有折叠自由度；不能仅依据0个可旋转键判为完全刚性。"),
    "F2_S27": (12, "羰基与亚砜氧竞争结合，中心环折叠需要检查。"),
    "F2_S28": (32, "11个可旋转键、多个NH/OH及羧酸；采用高采样预算，但须检查两性离子/质子化状态。"),
    "F3_N01": (8, "刚性单环N-氧化物，主要采样N-氧端的局部水/溶剂排列。"),
    "F3_N06": (8, "刚性稠合双环、一个主要N-氧受体区域。"),
    "F3_N14": (12, "N-氧与另一个环氮提供竞争受体区域。"),
    "F3_N15": (8, "溴代单环基本刚性；结构采样不由Br混合J基组设置替代。"),
    "F3_N16": (16, "羟甲基可向N-氧形成内氢键，也可被DMSO/水竞争打断。"),
    "F3_N17": (12, "稠环虽刚性，但N-氧与另一环氮的溶剂模式不同。"),
    "F3_N18": (12, "N-氧与含硫杂环的不同溶剂接近方向需要覆盖。"),
    "F3_N25": (12, "腈端和N-氧端的溶剂定位不同。"),
    "F3_N26": (12, "哌啶环椅式及N-氧取向不能用0个可旋转键概括。"),
    "F3_N28": (12, "环己基椅式及环外胺氧取向需要联合采样。"),
    "F3_N09": (12, "联吡啶环间扭转和两个不同受体位点。"),
    "F3_N10": (12, "刚性大稠环但存在N-氧/未氧化N两个结合区域。"),
    "F3_N19": (8, "刚性稠环、一个主要N-氧受体；重点为局部溶剂取向。"),
    "F3_N29": (24, "两个芳环间扭转加三个N-氧中心，溶剂排列组合较多。"),
    "F3_N30": (16, "羧酸、N-氧及环氮竞争氢键；酸/碱及互变异构状态需实验核实。"),
    "F3_N32": (16, "多环生物碱虽无普通可旋转键，仍需考虑饱和环形变及羰基/N-氧的溶剂分配。"),
    "F3_N33": (24, "柔性生物碱、醇、杂芳环与N-氧共存，具有多种竞争氢键/折叠。"),
    "F3_B01": (8, "小五元内酯；保留环起伏和两个氧周围溶剂取向。"),
    "F3_B03": (16, "七元内酯存在多种环构象；0个RDKit可旋转键不代表刚性。"),
    "F3_B05": (16, "取代七元内酯的环翻转、位阻和溶剂取向均需覆盖。"),
    "F3_B06": (12, "笼状内酯总体刚性，但扩环桥连骨架与两个氧的溶剂模式需检查。"),
    "F3_B08": (12, "芳氧酯/甲氧基扭转及三个氧受体位置不同。"),
    "F3_B14": (16, "七元环和苯基相对取向增加构象数。"),
    "F3_B15": (16, "多取代七元内酯的环形变及异丙基取向需采样。"),
    "F3_B16": (16, "螺环、七元内酯及四个氧提供多种环/溶剂组合。"),
    "F3_B18": (12, "桥环相对受限，但内酯结合位点和桥环形变需要复核。"),
    "F3_B20": (8, "芳环稠合五元内酯较刚性，主要是局部溶剂排列。"),
    "F3_B21": (16, "稠合环加OTBS的体积/扭转，需兼顾内酯端和硅醚端的溶剂位置。"),
    "F4_D01": (24, "5个可旋转键、二苯甲基与酰胺/亚砜竞争氢键；药物高覆盖采样。"),
    "F4_D02": (24, "丙基、氨基甲酸酯及两个NH；需覆盖侧链与多重氢键布局，并检查互变异构体。"),
    "F4_D03": (32, "两个杂环、多个甲氧基、NH及亚砜；多受体与柔性连接链采用增强预算。"),
    "F4_D04": (32, "杂环、柔性含氟醚链及NH/亚砜竞争氢键；同时需充分的含氟耦合构象覆盖。"),
    "F4_D05": (24, "羧酸、亚砜与芳基扭转组合；保留指定烯键几何，不混合不同几何异构体。"),
    "F4_D07": (24, "三环柔性侧链N-氧化物；主模型暂为净中性形式，真实盐型/质子化与溶剂必须核实。"),
}

OPTIONAL_CONTROLS = {
    "F2_E01": "简单环氧化物基线",
    "F2_E12": "OH与环氧氧竞争氢键",
    "F2_S01": "亚砜基线",
    "F2_S28": "高柔性/多氢键与潜在两性离子敏感性",
    "F3_N01": "N-氧化物强极性位点基线",
    "F3_B03": "中环内酯构象敏感性",
    "F4_D04": "复杂含氟药物与多位点氢键",
    "F4_D07": "复杂胺氧化物样品状态敏感性；先核实盐型",
}

ADVISORY_BANDS = {
    8: ("rigid_local_solvation", 12, 30),
    12: ("moderate_motif_or_ring", 20, 50),
    16: ("flexible_or_medium_ring", 30, 70),
    24: ("high_flexibility_or_drug", 50, 110),
    32: ("very_flexible_or_many_motifs", 70, 150),
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def plan(root: Path) -> tuple[dict, dict]:
    product_file = root / "registry/products.json"
    species_file = root / "registry/species.json"
    products = json.loads(product_file.read_text(encoding="utf-8"))
    species = json.loads(species_file.read_text(encoding="utf-8"))
    names = {p["compound_id"]: p["name"] for p in products}
    if set(names) != set(DECISIONS) or len(products) != 64:
        raise ValueError("The manually reviewed 64-product structure set changed: review decisions first")
    if len(species) != 67 or {s["compound_id"] for s in species} != set(names) | {"TMS"}:
        raise ValueError("Expected 66 product models and one TMS model; review changed model set")
    rows = []
    for s in species:
        cid = s["compound_id"]
        mol = Chem.MolFromSmiles(s["smiles"])
        if mol is None:
            raise ValueError(f"Invalid structure: {s['species_id']}")
        is_tms = cid == "TMS"
        n, reason = (8, "TMS无强氢键位点，但必须以相同溶剂簇组成独立采样并独立加权；不用产物权重代替。") if is_tms else DECISIONS[cid]
        tier, lo, hi = ADVISORY_BANDS[n]
        descriptors = {
            "rotatable_bonds_strict": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "hbond_donors": Descriptors.NumHDonors(mol),
            "hbond_acceptors": Descriptors.NumHAcceptors(mol),
            "heavy_atoms": mol.GetNumHeavyAtoms(),
            "rings": rdMolDescriptors.CalcNumRings(mol),
            "maximum_ring_size": max((len(x) for x in mol.GetRingInfo().AtomRings()), default=0),
            "formal_charge": Chem.GetFormalCharge(mol),
        }
        review_flags = []
        if not is_tms:
            review_flags += ["experimental_solvent_temperature_and_water_content_not_confirmed",
                             "computational_structure_not_experimental_identity_or_stereo_assignment"]
        if cid in {"F2_E11", "F2_E22"}:
            review_flags += ["separate_diastereomer_branches_do_not_pool_without_experimental_ratio"]
        if cid in {"F2_S28", "F3_N30", "F4_D02", "F4_D03", "F4_D04", "F4_D05", "F4_D07"}:
            review_flags += ["sample_state_protonation_or_tautomer_review_required"]
        if cid == "F4_D07":
            review_flags += ["neutral_DMSO_model_provisional_HCl_salt_and_CD3OD_branches_held"]
        jmethod = "NOT_APPLICABLE_REFERENCE" if is_tms else ("PBE_pcJ3" if cid.startswith("F4_") else "BENCHMARK_SELECT_PBE_pcJ3_OR_PBE0_pcJ2")
        rows.append({
            "species_id": s["species_id"], "compound_id": cid,
            "name": names.get(cid, "Tetramethylsilane (TMS reference)"),
            "smiles": s["smiles"], "model_role": s.get("role"),
            "scenario_id": SCENARIO, "main_solvent": "DMSO-d6",
            "explicit_DMSO_d6": 2, "explicit_H2O": 1, "explicit_D2O": 0,
            "outer_sampling_model": "ALPB(DMSO)", "outer_DFT_model": "CPCM(DMSO)",
            "temperature_K": 298.15, "temperature_status": "PROVISIONAL_UNTIL_EXPERIMENT_LOCKED",
            "solvent_composition_status": "MODELLING_ASSUMPTION_NOT_MEASURED_WATER_CONCENTRATION",
            "n_starts": n, "n_starts_per_species": n, "tier": tier,
            "suggested_unique_min": lo, "suggested_unique_max": hi,
            "unique_count_meaning": "ADVISORY_PLANNING_BAND_FOR_DISTINCT_RELEVANT_MINIMA_NOT_REQUIRED_OUTPUT_OR_TRUNCATION",
            "advisory_no_quota": True, "n_starts_is_not_n_final_conformers": True,
            "sampling_converged": False,
            "count_decision": "STRUCTURE_INFORMED_INITIAL_BUDGET_REVIEW_AFTER_FIRST_WAVE",
            "descriptors": descriptors, "reason": reason, "notes": reason,
            "review_flags": review_flags, "J_method": jmethod,
            "J_method_status": "USER_SELECTED_REQUIRES_VALIDATION_NOT_PROVEN_ACCURACY" if cid.startswith("F4_") else "REFERENCE_OR_BENCHMARK_DEPENDENT",
            "optional_control_recommended": cid in OPTIONAL_CONTROLS,
            "optional_control_reason": OPTIONAL_CONTROLS.get(cid),
            "optional_control_scenarios": ["DMSO_00_W00_continuum", "DMSO_02_W00_primary"] if cid in OPTIONAL_CONTROLS else [],
            "optional_controls_enabled": False,
        })
    main_rows = [r for r in rows if r["compound_id"] != "TMS"]
    aggregated = []
    for product in products:
        cid = product["compound_id"]
        branch_rows = [r for r in main_rows if r["compound_id"] == cid]
        first = branch_rows[0]
        aggregated.append({
            "compound_id": cid, "name": product["name"], "figure": int(cid[1]),
            "reaction_class": product["reaction_class"],
            "species_ids": [r["species_id"] for r in branch_rows],
            "n_species": len(branch_rows), "scenario_id": SCENARIO,
            "main_solvent": "DMSO-d6", "explicit_DMSO_d6": 2,
            "explicit_H2O": 1, "explicit_D2O": 0,
            "outer_sampling_model": "ALPB(DMSO)", "outer_DFT_model": "CPCM(DMSO)",
            "temperature_K": 298.15,
            "n_starts_per_species": {r["species_id"]: r["n_starts"] for r in branch_rows},
            "total_starts": sum(r["n_starts"] for r in branch_rows),
            "suggested_unique_min_per_species": first["suggested_unique_min"],
            "suggested_unique_max_per_species": first["suggested_unique_max"],
            "suggested_unique_min": first["suggested_unique_min"],
            "suggested_unique_max": first["suggested_unique_max"],
            "advisory_no_quota": True, "tier": first["tier"],
            "J_method": first["J_method"], "notes": first["reason"],
            "review_flags": first["review_flags"],
            "optional_control_recommended": first["optional_control_recommended"],
            "optional_controls_enabled": False,
        })
    summary = {
        "n_products": 64, "n_product_species": len(main_rows), "n_TMS_species": 1,
        "product_starts": sum(r["n_starts"] for r in main_rows),
        "TMS_starts": sum(r["n_starts"] for r in rows if r["compound_id"] == "TMS"),
        "total_starts": sum(r["n_starts"] for r in rows),
        "by_figure": {str(f): {"n_products": sum(r["figure"] == f for r in aggregated),
                               "n_starts": sum(r["total_starts"] for r in aggregated if r["figure"] == f)} for f in (2, 3, 4)},
        "product_species_by_start_count": {str(k): v for k, v in sorted(Counter(r["n_starts"] for r in main_rows).items())},
        "n_optional_control_products": len(OPTIONAL_CONTROLS),
        "n_enabled_control_starts": 0,
    }
    shared = {"schema_version": "4.0", "plan_version": PLAN_VERSION,
              "scenario_id": SCENARIO, "package_name": "工作2-计算资源包-（4）",
              "source_sha256": {"registry/products.json": sha(product_file), "registry/species.json": sha(species_file)},
              "summary": summary,
              "interpretation": "Initial independent GOAT start budgets only. Each start may produce many minima. Review rediscovery and final prediction convergence; do not impose the advisory minima band as a quota or cap.",
              "water_note": "One ordinary explicit H2O per cluster is a sensitivity/model assumption, not a measured bulk water content. DMSO continuum remains active concurrently.",
              "status": "PLANNED_NOT_QUANTUM_OPTIMIZED_NOT_CONVERGENCE_VALIDATED"}
    sampling = {**shared, "species": rows, "products": aggregated,
                "tms_reference": next(r for r in rows if r["compound_id"] == "TMS")}
    solvent = {**shared, "products": aggregated, "tms_reference": sampling["tms_reference"]}
    return sampling, solvent


def write_plan(root: Path) -> dict:
    sampling, solvent = plan(root)
    dump(root / "registry/sampling_plan_v4.json", sampling)
    dump(root / "registry/product_solvent_sampling_v4.json", solvent)
    s = sampling["summary"]
    text = f"""工作2-计算资源包-（4）：显式溶剂与GOAT数量依据

1. 主计算环境
全部64个产物均暂定：2个显式DMSO-d6＋1个普通H2O。
GOAT同时使用外围ALPB(DMSO)，DFT同时使用外围CPCM(DMSO)。
温度298.15 K为计划值。并未测量这些产物的实际样品溶剂、含水量、温度。
一个显式水不代表确定的体相含水量，不等于加入D2O；没有承诺所有交换性氢的实验化学位移准确。
不需要先运行纯连续介质或2 DMSO无水分支；两者是可选对照，不是前置步骤。

2. 数量到底是什么
本计划的n_starts是“不同初始溶质构象/显式溶剂布局的独立GOAT任务数”。
不是“最终只保留多少个构型”，也不是每个GOAT任务只输出一个构型。
每个GOAT任务可生成多个最低点；重复的最低点应去重，保留相关独特最低点进入后续流程。
suggested_unique_min/max只是评估计算资源的候选最低点参考带：不能凑数，也不能达到上限就截断。
真正的收敛应根据不同起点是否反复发现同一批低能构型，以及继续扩充后NMR预测是否稳定来判断。
若未覆盖关键环构象、氢键模式或发现新低能构型，应追加起点；若已稳定且独特最低点很少，不必硬凑参考带。
本版没有把达到起点数或最低点数自动视为科学收敛。

3. 人工决定的初始预算
8：较刚性、少量受体；12：多个局部氢键模式或受限环；16：柔性侧链/中环；
24：高柔性、多结合位点或复杂药物；32：长链、多供受体或更复杂药物。
使用RDKit可旋转键、环、供受体数量辅助，但不是按单一描述符机械分配。
尤其环辛烷、七元内酯、哌啶等，即使可旋转键为0，也需要专门环构象采样。
逐分子的理由与描述符见registry/sampling_plan_v4.json。

共64个产物，66个产品结构分支，另有1个TMS参照。
图2：{s['by_figure']['2']['n_products']}个产物，{s['by_figure']['2']['n_starts']}个独立起点。
图3：{s['by_figure']['3']['n_products']}个产物，{s['by_figure']['3']['n_starts']}个独立起点。
图4：{s['by_figure']['4']['n_products']}个产物，{s['by_figure']['4']['n_starts']}个独立起点。
产品总计{s['product_starts']}个，TMS另{s['TMS_starts']}个，总计{s['total_starts']}个主分支起点。
E11的endo/exo各8个，共16个；E22两个面异构模型各12个，共24个。
这些模型分别独立计算、独立加权；未经实验比例支持，不得混合为同一NMR样品预测。

4. 仍需人工核实的科学信息
所有结构是明确可计算的模型，并不等于已实验证实产物身份、绝对构型、溶剂或盐型。
F4_D07暂为净中性N-氧化物/DMSO模型；HCl盐、质子化位置及CD3OD分支仍不自动启用。
F2_S28的潜在两性离子状态、F3_N30的酸碱状态，以及部分苯并咪唑药物的互变异构需对照真实样品。
加水可能引起氢键重排甚至质子转移；不能把改变了化学状态的簇静默混入原结构权重。

5. J耦合
图4六个药物模型指定PBE/pcJ-3（内部标识PBE_pcJ3）。
图2/3保留PBE/pcJ-3与PBE0/pcJ-2的基准候选，未在本计划中擅自决定实验最佳方法。
指定方法不等于通过实验精度验证；TMS只承担同溶剂组成的化学位移参照，不要求J耦合。

6. 可选对照（本计划不自动生成或提交）
推荐8个代表物：F2_E01、F2_E12、F2_S01、F2_S28、F3_N01、F3_B03、F4_D04、F4_D07。
在样品状态明确后，可对照纯连续介质与2 DMSO无水模型。
如开展对照，每一种溶剂簇组成也必须有对应TMS参照，不能沿用2＋1组成的参照或权重。
这8个代表物是模型敏感性对照，不替代独立实验基准，也不自动证明其余产物的溶剂模型正确。
"""
    (root / "docs/构象数量依据.txt").write_text(text, encoding="utf-8")
    return s


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    print(json.dumps(write_plan(args.root.resolve()), ensure_ascii=False, indent=2))
