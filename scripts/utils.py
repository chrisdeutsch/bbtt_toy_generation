masspoints = [
    251, 260, 280, 300, 325, 350, 375, 400,
    450, 500, 550, 600, 700, 800, 900, 1000,
    1100, 1200, 1400, 1600
]


def addHeavyFlavourSplit(df):
    mask_zjets = (df["sample"] == "Z") | (df["sample"] == "Ztt")

    # B-jet truth matching
    b0, b1 = df["b0_truth_match"], df["b1_truth_match"]

    # Map taus to light
    df.loc[mask_zjets & (b0 == 15), "b0_truth_match"] = 0
    df.loc[mask_zjets & (b1 == 15), "b1_truth_match"] = 0

    mask_bb = (b0 == 5) & (b1 == 5)
    mask_bc = ((b0 == 5) & (b1 == 4)) | ((b0 == 4) & (b1 == 5))
    mask_cc = (b0 == 4) & (b1 == 4)
    mask_bl = ((b0 == 5) & (b1 == 0)) | ((b0 == 0) & (b1 == 5))
    mask_cl = ((b0 == 4) & (b1 == 0)) | ((b0 == 0) & (b1 == 4))
    mask_l = (b0 == 0) & (b1 == 0)

    df.loc[mask_zjets & mask_bb, "sample"] += "bb"
    df.loc[mask_zjets & mask_bc, "sample"] += "bc"
    df.loc[mask_zjets & mask_cc, "sample"] += "cc"
    df.loc[mask_zjets & mask_bl, "sample"] += "bl"
    df.loc[mask_zjets & mask_cl, "sample"] += "cl"
    df.loc[mask_zjets & mask_l, "sample"] += "l"
