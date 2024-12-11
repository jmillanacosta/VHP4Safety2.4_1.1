import pandas as pd

"""
This module provides functions to process gamma counter assay data for various substrates, cell types, and perturbagens (EDCs).

### Required Inputs:
- `raw_df`: A DataFrame containing the raw gamma counter data with columns such as `I-125 CPM`, `perturbagen`, `cell_type`, and `substrate`.
- `metadata`: A DataFrame specifying combinations of `cell_type` and `substrate_name` to process.
- `original_vol` and `new_vol`: Original and scaled assay volumes in microliters (floats).
- `edcs_dfs`: A dictionary mapping substrate names (keys) to their corresponding EDC response DataFrames.
- `bg`: A DataFrame containing background responses for substrates.

### Outputs:
- Consolidated DataFrame combining results for all combinations, with columns for `cell_type` and `substrate_name`.
"""


def tag_rows(raw_df):
    """
    Placeholder function to annotate raw_df with perturbagen, cell type, and substrate.
    This function should be implemented to add relevant annotations to the input DataFrame.

    :param raw_df: pd.DataFrame: The raw input DataFrame.
    :return: pd.DataFrame: Annotated DataFrame.
    """
    raise NotImplementedError(
        "This function is a placeholder and should be implemented."
    )


def get_edcs_cpms(raw_df, substrate, cell_type, original_vol, new_vol):
    """
    Calculate CPMs for EDCs for a specific substrate and cell type.

    :param raw_df: pd.DataFrame: Input DataFrame.
    :param substrate: str: The substrate to analyze.
    :param cell_type: str: The cell type to analyze.
    :param original_vol: float: Original assay volume in microliters.
    :param new_vol: float: New volume for scaling CPM values, in microliters.
    :return: pd.DataFrame: A DataFrame with CPM values for each perturbagen.
    """
    volume_ratio = new_vol / original_vol
    filtered_df = raw_df[
        (raw_df["substrate"] == substrate) & (raw_df["cell_type"] == cell_type)
    ]

    cpms = {
        perturbagen: {"cpm1": "", "cpm2": "", "cpm_average": "", f"cpm_{new_vol}": ""}
        for perturbagen in filtered_df["perturbagen"].unique()
    }

    for perturbagen in cpms:
        perturbagen_cpms = filtered_df[filtered_df["perturbagen"] == perturbagen][
            "I-125 CPM"
        ].values
        if len(perturbagen_cpms) >= 2:
            cpms[perturbagen]["cpm1"] = perturbagen_cpms[0]
            cpms[perturbagen]["cpm2"] = perturbagen_cpms[1]
            cpms[perturbagen]["cpm_average"] = (
                perturbagen_cpms[0] + perturbagen_cpms[1]
            ) / 2
            cpms[perturbagen][f"cpm_{new_vol}"] = (
                cpms[perturbagen]["cpm_average"] * volume_ratio
            )

    return pd.DataFrame(cpms).transpose()


def get_substrate_cpms(raw_df, substrate_names, original_vol, new_vol):
    """
    Calculate CPMs for specific substrates.

    :param raw_df: pd.DataFrame: Input DataFrame.
    :param substrate_names: list: List of substrate names to analyze.
    :param original_vol: float: Original assay volume in microliters.
    :param new_vol: float: New volume for scaling CPM values, in microliters.
    :return: pd.DataFrame: A DataFrame with CPM values for each substrate.
    """
    volume_ratio = new_vol / original_vol
    filtered_df = raw_df[raw_df["perturbagen"].isin(substrate_names)]

    cpms = {
        perturbagen: {"cpm1": "", "cpm2": "", "cpm_average": "", f"cpm_{new_vol}": ""}
        for perturbagen in filtered_df["perturbagen"].unique()
    }

    for perturbagen in cpms:
        perturbagen_cpms = filtered_df[filtered_df["perturbagen"] == perturbagen][
            "I-125 CPM"
        ].values
        if len(perturbagen_cpms) >= 2:
            cpms[perturbagen]["cpm1"] = perturbagen_cpms[0]
            cpms[perturbagen]["cpm2"] = perturbagen_cpms[1]
            cpms[perturbagen]["cpm_average"] = (
                perturbagen_cpms[0] + perturbagen_cpms[1]
            ) / 2
            cpms[perturbagen][f"cpm_{new_vol}"] = (
                cpms[perturbagen]["cpm_average"] * volume_ratio
            )

    return pd.DataFrame(cpms).transpose()


def get_subs_edcs_cell_cpms(  #TODO figure out role of DMSO as bg
    raw_df, substrate_name, original_vol, new_vol, cell_type, edcs_df, substrate_df
):
    """
    Calculate corrected CPMs for a specific substrate, cell type, and EDC combinations.

    :param raw_df: pd.DataFrame: Input DataFrame.
    :param substrate_name: str: Substrate to analyze.
    :param original_vol: float: Original volume of the assay in microliters.
    :param new_vol: float: New volume for scaling CPMs.
    :param cell_type: str: Cell type to analyze.
    :param edcs_df: pd.DataFrame: EDCs response DataFrame.
    :param substrate_df: pd.DataFrame: Substrate background response DataFrame.
    :return: pd.DataFrame: A DataFrame with corrected CPMs.
    """
    df = edcs_df[[f"cpm_{new_vol}"]].copy()
    df.columns = ["tracermix"]
    filtered_raw_df = raw_df[
        (raw_df["substrate"] == substrate_name) & (raw_df["cell_type"] == cell_type)
    ]

    df = df.join(
        get_edcs_cpms(
            raw_df,
            substrate_name,
            cell_type,
            original_vol,
            new_vol,
        ),
        lsuffix="_left",
        rsuffix="_right",
    )
    df = df.replace("", None).astype(float)
    df["cpm_1_corrected"] = df["cpm1"] / df["tracermix"]
    df["cpm_2_corrected"] = df["cpm2"] / df["tracermix"]
    df["average_cpm_corrected"] = df[["cpm_1_corrected", "cpm_2_corrected"]].mean(
        axis=1
    )
    df["uptake 1 tov DMSO"] = (
        df["cpm_1_corrected"] / df.loc["DMSO 0.5%", "average_cpm_corrected"]
    )
    df["uptake 2 tov DMSO"] = (
        df["cpm_2_corrected"] / df.loc["DMSO 1%", "average_cpm_corrected"]
    )
    df["avg uptake"] = df[["uptake 1 tov DMSO", "uptake 2 tov DMSO"]].mean(axis=1)
    df["std uptake"] = df[["uptake 1 tov DMSO", "uptake 2 tov DMSO"]].std(axis=1)
    return df


def run_all_combinations(raw_df, metadata, original_vol, new_vol, edcs_dfs, bg):
    """
    Automates the processing of CPM calculations for all cell/substrate combinations
    based on metadata and returns a consolidated DataFrame.

    :param raw_df: pd.DataFrame: Input DataFrame with raw gamma counter assay data.
    :param metadata: pd.DataFrame: Metadata file with columns "cell_type" and "substrate_name".
    :param original_vol: float: Original assay volume in microliters.
    :param new_vol: float: New volume for scaling CPM values, in microliters.
    :param edcs_dfs: dict: A dictionary mapping substrate names to their corresponding EDC response DataFrames.
    :param bg: pd.DataFrame: Background responses for substrates.
    :return: pd.DataFrame: A consolidated DataFrame with results for all combinations,
        including `cell_type` and `substrate_name` columns.
    """
    combined_results = []

    for _, row in metadata.iterrows():
        cell_type = row["cell_type"]
        substrate_name = row["substrate_name"]
        edcs_df = edcs_dfs.get(substrate_name)

        if edcs_df is None:
            raise ValueError(
                f"No EDC response data found for substrate: {substrate_name}"
            )

        result_df = get_subs_edcs_cell_cpms(
            raw_df, substrate_name, original_vol, new_vol, cell_type, edcs_df, bg
        )

        result_df["cell_type"] = cell_type
        result_df["substrate_name"] = substrate_name

        combined_results.append(result_df)

    return pd.concat(combined_results, ignore_index=True)
