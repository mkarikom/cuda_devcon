import os
import scanpy as sc
import pandas
import pickle
# Get the data:
# 1. use GEOparse to get the supplementary file FTP handles
# 2. unzip and read this data into MultiVI

def load_inflam_new(datadir,metadir,serializedir):
    phenotable = pandas.read_csv(f"{metadir}/phenotable.csv")
    phenotable.index = phenotable['shortname']

    if not os.path.exists(f"{serializedir}/gact_dict.pickle"):
        if os.path.exists(f"{datadir}"):
            adata_dict = dict()
            gact_dict = dict()
            for pid in phenotable.index:
                print(pid)
                # read in the multiome data
                if not os.path.isdir(f"{unzipdatadir}/{pid}"):
                    os.makedirs(f"{unzipdatadir}/{pid}")
                    for fn in os.listdir(f"{datadir}/{pid}"):
                        print(f"unzipping {fn} from {pid}")
                        newname = os.path.splitext(fn)[0]
                        os.system(f"gunzip -ck {datadir}/{pid}/{fn} > {unzipdatadir}/{pid}/{newname}")
                else:
                    print(f"{unzipdatadir}/{pid} already exists, using previously extracted data")
                print(f"importing signac-processed multiome data from {pid}")
                adata = scvi.data.read_10x_multiome(f"{unzipdatadir}/{pid}")
                # make sure the barcodes are unique when eventually concatenated
                adata.obs['barcode'] = [i+"-"+pid for i in adata.obs.index.tolist()]
                adata.obs = adata.obs.set_index('barcode')
                # add the phenodata to the obs table
                adata.var_names_make_unique()
                adata.obs['disease'] = phenotable.loc[pid].disease
                adata.obs['ifn1'] = phenotable.loc[pid].ifn1
                adata.obs['pid'] = pid
                # read in the predicted gene activity
                gact = pandas.read_csv(f"{metadir}/{pid}/geneactivity.csv")
                # set the index of GACT to the barcode 
                gact = gact.set_index('Unnamed: 0',drop=True)
                # get rid of any trailing decimels in the barcode names
                gact.columns = gact.columns.str.replace(".1","")
                # # only keep GACT data that overlaps with GEX
                # gact = gact.loc[gact.index.intersection(adata.var.index),gact.columns.intersection(adata.obs.index.array)]
                # # only keep multiome cells for which GACT data is available
                # adata = adata[gact.columns.intersection(adata.obs.index.array)]
                adata_dict[pid] = adata
                gact = sc.AnnData(gact.transpose(),obs=gact.columns.to_frame(name="barcode"),var=gact.index.to_frame(name="gid"))
                gact.var_names_make_unique()
                # make sure the barcodes are unique when eventually concatenated
                gact.obs['barcode'] = [i+"-"+pid for i in gact.obs.index.tolist()]
                gact.obs = gact.obs.set_index('barcode')
                # add the phenodata to the obs table
                gact.obs['disease'] = phenotable.loc[pid].disease
                gact.obs['ifn1'] = phenotable.loc[pid].ifn1
                gact.obs['pid'] = pid
                gact_dict[pid] = gact

            if not os.path.exists(serializedir):
                os.makedirs(serializedir)
            with open(f"{serializedir}/gact_dict.pickle", 'wb') as handle:
                pickle.dump(gact_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
            with open(f"{serializedir}/adata_dict.pickle", 'wb') as handle:
                pickle.dump(adata_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            print("multiome or phenodata missing")
    else:
        print(f"loading serialized gex from {serializedir}/adata_dict.pickle")
        with open(f"{serializedir}/adata_dict.pickle", 'rb') as handle:
            adata_dict = pickle.load(handle)

        print(f"loading serialized gact from {serializedir}/gact_dict.pickle")
        with open(f"{serializedir}/gact_dict.pickle", 'rb') as handle:
            gact_dict = pickle.load(handle)
    
    return dict(adata=adata_dict,gact=gact_dict)