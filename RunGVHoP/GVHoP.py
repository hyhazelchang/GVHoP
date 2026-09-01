#!/usr/bin/python3

# GVHoP.py

# Hsin-Ying Chang <hyhazelchang@gmail.com>
# v1 2026/08/28
# from load predictor

# cd /home/xinchang/projects/girush02/girush02.26/DATA/pyscripts/
# python3 GVHoP.py --GVOGs_in=../source_data/example_inputs/ex_GVOGs.tsv --GVEUKs_in=../source_data/example_inputs/ex_GVEUKs.tsv --sample_ls=../source_data/example_inputs/run1_sample.ls --out_dir=../outdir/

import os
import argparse
import numpy as np
import pandas as pd
from functools import reduce
from collections import defaultdict
import joblib
import xgboost as xgb
from sklearn.preprocessing import StandardScaler
import torch
from torch.utils.data import Dataset
import torch.nn as nn
from torch.utils.data import DataLoader
from sklearn.compose import ColumnTransformer
from torch.utils.data import Sampler, RandomSampler, SequentialSampler
from GVHoP_NN import MetaDataset
from GVHoP_NN import MetaNN
from GVHoP_NN import SafeBatchSampler

def run_XGBclf(data_df, scaler_col, model, prob_out, hosts):
    # Load new data (unlabeled)
    X_test = data_df.to_numpy()
    X_test_final = []
    for i, scaler in enumerate(model.named_steps['scalers']):
        X_block = X_test[:, i*scaler_col:(i+1)*scaler_col]
        X_block_scaled = scaler.transform(X_block)
        X_test_final.append(X_block_scaled)
    X_test_scaled = np.concatenate(X_test_final, axis=1)
    X_test_selected = model.named_steps['feature_selector'].transform(X_test_scaled)    
    
    # Predict
    # y_pred = model.named_steps['classifier'].predict(X_test_selected)
    # y_predpr = model.named_steps['classifier'].predict_proba(X_test_selected)
    y_scores = model.named_steps['classifier'].predict(X_test_selected, output_margin=True)

    # Output the prediction results
    # prob_out['pred_label'] = y_pred
    prob_out_new = pd.concat([prob_out.reset_index(drop=True), 
                              pd.DataFrame(y_scores, columns = hosts).reset_index(drop=True)], 
                              axis=1)
    return prob_out_new

def aggregate_prob(probs, node_dict, node_name):
    return sum(probs[i] for i in node_dict[node_name])

def climbing_inference(prob_tb, hierarch_list, node_dict, threshold):
    """
    Function for climbing inference for HSC
    # ex: test_label: 2130 (from most specific → most general)
    # pred_label: --30 (from most specific → most general)
    # from predicted leaf to root
    """
    # ----------------------------
    # Reverse the hierarch from most specific → most general
    # ----------------------------
    rev_hierarch_list = [sublist[::-1] for sublist in hierarch_list]

    # ----------------------------
    # Define correct answers
    # ----------------------------
    prob_out_levels = {}
    prob_out_levels["testset"] = prob_tb["testset"].to_list()
    # prob_out_levels["test_label"] = [[x] for x in prob_tb["test_label"].to_list()]
    y_pred = prob_tb["pred_label"].to_numpy()
    y_predpr = prob_tb.iloc[:, 2:].to_numpy()

    hierarch_dict = {}
    for i in range(len(hierarch_list[0])):
        labels = sorted({sublist[i] for sublist in rev_hierarch_list})
        hierarch_dict[f"level{i}"] = labels
    # ----------------------------
    # Climbing inference rule with threshold
    # ----------------------------
    prob_out_levels["pred_label"] = []
    prob_out_levels["pred_hierarchy"] = []
    for i in range(len(hierarch_list[0])):
        prob_out_levels[f"level{i}_prob"] = []

    for j in range(len(y_pred)):
        pred_l = y_pred[j] # predicted label: 0, 1, 2...
        predpr = y_predpr[j] # probabilities of predictions
        prediction = hierarch_dict["level0"][pred_l] # predicted name
        hierach_of_pred_l = [item for sublist in rev_hierarch_list if prediction in sublist for item in sublist] # hierarch of predicted label
        pred_label_index = []
        pred_label_hierarch = []
        # climbing
        for h in range(len(hierach_of_pred_l)):
            p = aggregate_prob(predpr, node_dict, hierach_of_pred_l[h])
            if p >= threshold:
                prob_out_levels[f"level{h}_prob"].append(p)
                pr = hierarch_dict[f"level{h}"].index(hierach_of_pred_l[h])
                pred_label_index.append(str(pr))
                pred_label_hierarch.append(hierach_of_pred_l[h])
            else:
                prob_out_levels[f"level{h}_prob"].append("-")
                pred_label_index.append("-")
                pred_label_hierarch.append("-")
        prob_out_levels["pred_label"].append(",".join(pred_label_index))
        prob_out_levels["pred_hierarchy"].append(",".join(pred_label_hierarch))
    df = pd.DataFrame(prob_out_levels)
    return df

def count_predictions_BU(prob_tb, all_pred, samples_pred):
    """Count hierarchy outputs"""
    pred_roots = 0
    test_ids = prob_tb["testset"].to_list()
    pred_label = [h.split(",") for h in prob_tb["pred_label"].to_list()]
    pred_hierarch = [h.split(",") for h in prob_tb["pred_hierarchy"].to_list()]
    for i in range(len(pred_label)):
        id = test_ids[i]
        for j in range(len(pred_label[0])):
            if pred_hierarch[i][j] != "-":			
                for k in range(j, len(pred_label[0])):
                    samples_pred[id][k][pred_hierarch[i][k]] += 1
                all_pred[id][pred_hierarch[i][k]] += 1
                break
            elif j == len(pred_label[0]) - 1:
                all_pred[id]["root"] += 1
                pred_roots += 1

def flatten_nested(lst):
    flattened_lst = []
    for item in lst:
        if isinstance(item, list):
            for subitem in flatten_nested(item):
                if subitem not in flattened_lst:
                    flattened_lst.append(subitem)
        else:
            if item not in flattened_lst:
                flattened_lst.append(item)
    return flattened_lst    

def main():
    parser = argparse.ArgumentParser(
            description=("Giant Virus-Host Predictor"),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--GVOGs_in",
                        type=str,
                        default=None,
                        help="The directory of input GVOGs files.")
    parser.add_argument("--GVEUKs_in",
                        type=str,
                        default=None,
                        help="The directory of input GVEUKs files.")
    parser.add_argument("--sample_ls",
                        type=str,
						required=True,
                        default=None)
    parser.add_argument("--out_dir",
                        type=str,
						required=True,
                        default='./outdir/',
                        help="The directory of output files.") 
       
    # Defining variables from input
    args = parser.parse_args()
    GVOGs_in = args.GVOGs_in
    GVEUKs_in = args.GVEUKs_in
    sample_ls = [line.strip() for line in open(args.sample_ls, "r")]
    out_dir = args.out_dir

    # Create output directory
    os.makedirs(out_dir, exist_ok=True)

    # Host lists
    top_ls = ['Algae', 'Amoeba', 'Fungi', 'Heteroflagellate', 'Metazoa']
    intermediate_ls = ['Amoeboid', 'Amoebozoa', 'Blastocladiomycota', 'Brownalgae', 'Chlorophyta', 'Chytridiomycota', 'Haptophyta', 'Invertebrate', 'Nonamoeboid', 'Otheralgae', 'Vertebrate']
    bottom_ls = ['Aves', 'Bathycoccaceae', 'Blastocladiomycetes', 'Chlorellaceae', 'Choanocafe', 'Chytridiomycetes', 'Coccolithophyceae', 'Discoba', 'Discosea', 'Insecta', 'Malacostraca', 'Mamiellaceae', 'Mammalia', 'Otherchlorophyta', 'Otherhaptophyta', 'Otherinvertebrate', 'Othervertebrate', 'Phaeophyceae', 'Prymnesiaceae', 'Stramenopile', 'Tubulinea']

    # Create the dictionary for saving the raw scores
    prob_out_dict = defaultdict(list)

    """ clf_GVOGs """
    # Parse data from input file
    data_df = pd.read_csv(GVOGs_in, index_col=0, header=0, sep='\t')
    data_df.index.name = "testset"
    # Select the samples
    data_df = data_df.reindex(sample_ls, fill_value=0)
    # Select GVHoP all features
    feature_cols = open("../source_data/features/GVHoP_GVOGs_all.tsv").readline().strip().split('\t')
    data_df = data_df.reindex(columns=feature_cols, fill_value=0)
    # Select GVHoP top features
    top_feature_cols = open("../source_data/features/GVHoP_GVOGs_top.tsv").readline().strip().split('\t')
    top_GVOGs_df = data_df.reindex(columns=top_feature_cols, fill_value=0)
    # Create a dataframe for saving probability
    samples = data_df.index.tolist()
    prob_out = pd.DataFrame({'testset': samples})
    # Load model
    print("Initialize clf_GVOGs...")
    for i in range(1, 101):
        ## clf_top ##
        print(f"clf_GVOGs_top_{i}_model")
        model = joblib.load(f"../XGBclf/clf_GVOGs/top/XGB_{i}.joblib")
        host_ls = [f"gc_s_{host}" for host in top_ls]
        prob_out_top = run_XGBclf(data_df, 8293, model, prob_out, host_ls)
        ## clf_intermediate ##
        print(f"clf_GVOGs_intermediate_{i}_model")
        model = joblib.load(f"../XGBclf/clf_GVOGs/intermediate/XGB_{i}.joblib")
        host_ls = [f"gc_s_{host}" for host in intermediate_ls]
        prob_out_intermediate = run_XGBclf(data_df, 8293, model, prob_out, host_ls)
        ## clf_bottom ##
        print(f"clf_GVOGs_bottom_{i}_model")
        model = joblib.load(f"../XGBclf/clf_GVOGs/bottom/XGB_{i}.joblib")
        host_ls = [f"gc_s_{host}" for host in bottom_ls]
        prob_out_bottom = run_XGBclf(data_df, 8293, model, prob_out, host_ls)
        #
        prob_out_dfs = [prob_out_top, prob_out_intermediate, prob_out_bottom]
        # merge the dataframes on the testset
        prob_out_new = reduce(lambda left, right: pd.merge(left, right, on='testset', how='outer'), prob_out_dfs)
        prob_out_dict[i].append(prob_out_new)
    print("Finish clf_GVOGs!")

    """ clf_GVEUKs """
    # Parse data from input file
    data_df = pd.read_csv(GVEUKs_in, index_col=0, header=0, sep='\t')
    data_df.index.name = "testset"
    # Select the samples
    data_df = data_df.reindex(sample_ls, fill_value=0)
    # Select GVHoP all features
    feature_cols = open("../source_data/features/GVHoP_GVEUKs_all.tsv").readline().strip().split('\t')
    data_df = data_df.reindex(columns=feature_cols, fill_value=0)
    # Select GVHoP top features
    top_feature_cols = open("../source_data/features/GVHoP_GVEUKs_top.tsv").readline().strip().split('\t')
    top_GVEUKs_df = data_df.reindex(columns=top_feature_cols, fill_value=0)
    # Create a dataframe for saving probability
    samples = data_df.index.tolist()
    prob_out = pd.DataFrame({'testset': samples})
    # Load model
    print("Initialize clf_GVEUKs...")
    for i in range(1, 101):
        ## clf_top ##
        print(f"clf_GVEUKs_top_{i}_model")
        model = joblib.load(f"../XGBclf/clf_GVEUKs/top/XGB_{i}.joblib")
        host_ls = [f"hgt_s_{host}" for host in top_ls]
        prob_out_top = run_XGBclf(data_df, 57250, model, prob_out, host_ls)
        ## clf_intermediate ##
        print(f"clf_GVEUKs_intermediate_{i}_model")
        model = joblib.load(f"../XGBclf/clf_GVEUKs/intermediate/XGB_{i}.joblib")
        host_ls = [f"hgt_s_{host}" for host in intermediate_ls]
        prob_out_intermediate = run_XGBclf(data_df, 57250, model, prob_out, host_ls)
        ## clf_bottom ##
        print(f"clf_GVEUKs_bottom_{i}_model")
        model = joblib.load(f"../XGBclf/clf_GVEUKs/bottom/XGB_{i}.joblib")
        host_ls = [f"hgt_s_{host}" for host in bottom_ls]
        prob_out_bottom = run_XGBclf(data_df, 57250, model, prob_out, host_ls)
        #
        prob_out_dfs = [prob_out_top, prob_out_intermediate, prob_out_bottom]
        # merge the dataframes on the testset
        prob_out_new = reduce(lambda left, right: pd.merge(left, right, on='testset', how='outer'), prob_out_dfs)
        prob_out_dict[i].append(prob_out_new)
    print("Finish clf_GVEUKs!")

    """ GVHoP_clf_integration """
    # Integrate the feature sets
    top_feat = top_GVOGs_df.merge(top_GVEUKs_df, on='testset', how='left')
    for sample, dfs in prob_out_dict.items():
        meta_df_RawScores = dfs[0].merge(dfs[1], on="testset", how="left")
        meta_df_feat = meta_df_RawScores.merge(top_feat, on="testset", how="left")
        prob_out_dict[sample] = meta_df_feat
    meta_df = pd.concat(prob_out_dict.values(), ignore_index=True)
    #
    # Get feature arrays
    FEATURES = meta_df.columns[1:]
    labels = bottom_ls
    n_classes_global = 21
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("Initialize meta-model... ")
    # 0. load GVMAGs testing set
    X_test = meta_df[FEATURES]
    #
    # 1. Load the preprocessor back into memory
    loaded_preprocessor = joblib.load("../NNclf/preprocessor_all.joblib")
    #
    # 2. Transform new data (CRITICAL: Use .transform(), NOT .fit_transform())
    X_test_scaled = loaded_preprocessor.transform(X_test)
    X_test_scaled[np.isnan(X_test_scaled)] = 0
    #
    # 3. make a df to store the data
    prob_out = pd.DataFrame(meta_df["testset"])
    #
    # 4. load dataset
    test_dataset = torch.tensor(X_test_scaled, dtype=torch.float32)
    safe_sampler = SafeBatchSampler(test_dataset, batch_size=8, shuffle=False)
    test_loader = DataLoader(test_dataset, batch_sampler=safe_sampler, num_workers=2, pin_memory=True)
    #
    # 5. initiate same model architecture
    model_path = "../NNclf/nn_all.pt"
    model = MetaNN(n_classes=n_classes_global)
    model.to(device)
    model.load_state_dict(torch.load(model_path, weights_only=True))
    #
    # 6. CRITICAL: Set the model to evaluation mode
    model.eval()
    probs = []
    preds = []
    with torch.no_grad():
        for xb in test_loader:
            xb = xb.to(device)
            logits = model(xb)
            # Apply Softmax to convert logits to probabilities
            probabilities = torch.softmax(logits, dim=1)
            probs.append(probabilities.cpu().numpy())
            predictions = torch.argmax(logits, dim=1)
            preds.append(predictions.cpu().numpy())
    prob_out["pred_label"] = np.concatenate(preds)
    prob_out[labels] = np.vstack(probs)
    prob_out.to_csv(f"{out_dir}prob_all.tsv", sep="\t", index=False, header=True)
    print("Finish meta-model! ")

    """ Get prediction hierarchy """
    hierarch_list = [['Metazoa', 'Vertebrate', 'Aves'], ['Algae', 'Chlorophyta', 'Chlorellaceae'], ['Amoeba', 'Amoebozoa', 'Discosea'], ['Heteroflagellate', 'Amoeboid', 'Discoba'], ['Algae', 'Chlorophyta', 'Otherchlorophyta'], ['Algae', 'Haptophyta', 'Coccolithophyceae'], ['Algae', 'Haptophyta', 'Otherhaptophyta'], ['Metazoa', 'Vertebrate', 'Mammalia'], ['Heteroflagellate', 'Nonamoeboid', 'Choanocafe'], ['Algae', 'Otheralgae', 'Stramenopile'], ['Algae', 'Chlorophyta', 'Bathycoccaceae'], ['Algae', 'Brownalgae', 'Phaeophyceae'], ['Fungi', 'Chytridiomycota', 'Chytridiomycetes'], ['Metazoa', 'Invertebrate', 'Insecta'], ['Metazoa', 'Invertebrate', 'Otherinvertebrate'], ['Amoeba', 'Amoebozoa', 'Tubulinea'], ['Algae', 'Chlorophyta', 'Mamiellaceae'], ['Metazoa', 'Vertebrate', 'Othervertebrate'], ['Fungi', 'Blastocladiomycota', 'Blastocladiomycetes'], ['Metazoa', 'Invertebrate', 'Malacostraca'], ['Algae', 'Haptophyta', 'Prymnesiaceae']]
    node_dict = {'Algae': [13, 3, 6, 19, 1, 14, 11, 18, 17], 'Chlorophyta': [13, 3, 1, 11], 'Otherchlorophyta': [13], 'Chlorellaceae': [3], 'Fungi': [5, 2], 'Chytridiomycota': [5], 'Chytridiomycetes': [5], 'Metazoa': [12, 10, 16, 15, 0, 9], 'Vertebrate': [12, 16, 0], 'Mammalia': [12], 'Invertebrate': [10, 15, 9], 'Malacostraca': [10], 'Othervertebrate': [16], 'Amoeba': [20, 8], 'Amoebozoa': [20, 8], 'Tubulinea': [20], 'Haptophyta': [6, 14, 18], 'Coccolithophyceae': [6], 'Heteroflagellate': [7, 4], 'Amoeboid': [7], 'Discoba': [7], 'Otheralgae': [19], 'Stramenopile': [19], 'Bathycoccaceae': [1], 'Otherhaptophyta': [14], 'Mamiellaceae': [11], 'Nonamoeboid': [4], 'Choanocafe': [4], 'Otherinvertebrate': [15], 'Blastocladiomycota': [2], 'Blastocladiomycetes': [2], 'Prymnesiaceae': [18], 'Discosea': [8], 'Aves': [0], 'Insecta': [9], 'Brownalgae': [17], 'Phaeophyceae': [17]}
    # bottom-up
    prob_out_levels_BU = climbing_inference(prob_out, hierarch_list, node_dict, 0.75)
    prob_out_levels_BU.to_csv(f"{out_dir}h_prob_all.tsv", sep="\t", index=False, header=True)

    # Create a dictionary for storing predictions
    all_pred_BU = defaultdict(dict)
    sorted_hierarch_list = sorted(hierarch_list)
    flat_hierarch_list = flatten_nested(sorted_hierarch_list)
    for label in flat_hierarch_list:
        for id in sample_ls:
            all_pred_BU[id][label] = 0
    ## add a root
    for id in sample_ls:
        all_pred_BU[id]["root"] = 0
    
    # Create a dictionary for storing predictions of listed samples
    samples_pred_BU = defaultdict(lambda: defaultdict(dict))
    rev_hierarch_list = [sublist[::-1] for sublist in sorted_hierarch_list]
    for i in range(len(hierarch_list[0])):
        labels = [sublist[i] for sublist in rev_hierarch_list]
        for id in set(prob_out["testset"]): ##sample_ids
            for label in labels:
                samples_pred_BU[id][i][label] = 0
		
    # Count predictions
    count_predictions_BU(prob_out_levels_BU, all_pred_BU, samples_pred_BU)

    """ Output counts by levels """
    by_level = defaultdict(dict)
    for sample_id in samples_pred_BU:
        for level in samples_pred_BU[sample_id]:
            by_level[level][sample_id] = samples_pred_BU[sample_id][level]
    # Save each k to a file
    for level, data in by_level.items():
        df = pd.DataFrame(data).T
        df.to_csv(f"{out_dir}pred_out_level{level}.tsv", sep='\t', index=True, header=True)

if __name__ == '__main__':
	main()