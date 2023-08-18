import pandas as pd
import json
from neo4j import GraphDatabase
import random
import time
import os
import glob
import math
from collections import Counter

csv_path = '/home/johannes/Desktop/db_test_data/csv/'
methylation_path = '/home/johannes/Desktop/db_test_data/methylation.csv'
mutation_path = '/home/johannes/Desktop/db_test_data/somatic_mutations_Unified.csv'
network_paths = ['/home/johannes/Desktop/db_test_data/THCA/dysregnet/tpm-exp.fea']
stats_paths = ['/home/johannes/Desktop/db_test_data/THCA/dysregnet/tpm-exp-stats.csv']
auth = ('neo4j', '12345678')
uri = 'bolt://localhost:7687'

def get_gene_ids(keys):
    gene_set = set()
    for key in keys:
        genes = key.split(', ')
        gene1 = genes[0].strip("()'")
        gene2 = genes[1].strip("()'")
        gene_set.add(gene1)
        gene_set.add(gene2)
    return gene_set


def write_patient_csv(patient_ids):
    with open(csv_path + "Patient.csv", "w") as f:
        f.write("patient_id,mapping_id\n")
        for key in patient_ids.keys():
            f.write(",".join([patient_ids[key], key]))
            f.write("\n")


def write_gene_csv(genes, methylation, mutation):

    with open(csv_path + "Gene.csv", "w") as f:
        f.write("gene_id,mutation,methylation\n")
        for gene in genes:
            gene_mutation = "%.5f" % mutation[gene] if gene in mutation else 0
            gene_methylation = "%.5f" % methylation[gene] if gene in methylation else ""
            f.write(",".join((gene, str(gene_mutation), str(gene_methylation))))
            f.write("\n")


def write_regulation_csv(map, edge_types):
    with open(csv_path + "Regulation.csv", "w") as f:
        f.write("regulation_id,source,target,fraction,direction\n")
        for key in map.keys():
            edge_type = edge_types[key] 
            genes = key.split(', ')
            source = genes[0].strip("()'")
            target = genes[1].strip("()'")
            regulation_id = source + ":" + target
            fraction = "%.5f" % (sum([value != 0 for value in map[key].values()]) / len(map[key].keys()))
            direction = "+" if edge_type > 0 else "-"
            f.write(",".join((regulation_id, source, target, str(fraction), direction)))
            f.write("\n")


def write_dysregulated_csv(map):
    with open(csv_path + "DYSREGULATED.csv", "w") as f:
        f.write("mapping_id,regulation_id,value\n")

        for key in map.keys():

            genes = key.split(', ')
            source = genes[0].strip("()'")
            target = genes[1].strip("()'")
            regulation_id = source + ":" + target
            values = map[key]

            for mapping_id in values.keys():
                if (values[mapping_id] != 0):
                    f.write(",".join((mapping_id, regulation_id, str(values[mapping_id]))))
                    f.write("\n")

def write_methylation_csv(methylation):
    with open(csv_path + "METHYLATED.csv", "w") as f:
        f.write("gene_id,patient_id,methylation\n")
        for rowIndex, row in methylation.iterrows():  # iterate over rows
            for columnIndex, value in row.items():
                f.write(",".join((rowIndex,columnIndex,str("%.5f" % value))))
                f.write("\n")

def write_mutation_csv(mutation_pairs):
    with open(csv_path + "MUTATED.csv", "w") as f:
        f.write("gene_id,patient_id\n")
        for gene, patient in mutation_pairs:
            f.write(",".join((gene,patient)))
            f.write("\n")

def get_pre_query_list():
    query_list = [
        f"CREATE CONSTRAINT FOR (p:Patient) REQUIRE p.patient_id IS UNIQUE;",
        f"CREATE CONSTRAINT FOR (c:Cancer) REQUIRE c.cancer_id IS UNIQUE;",
    ]
    return query_list

def get_pre_query_list_cancer(cancer):
    query_list = [
        f"CREATE CONSTRAINT FOR (p:{cancer}_Patient) REQUIRE p.patient_id IS UNIQUE;",
        f"CREATE CONSTRAINT FOR (p:{cancer}_Patient) REQUIRE p.mapping_id IS UNIQUE;",
        f"CREATE CONSTRAINT FOR (r:{cancer}_Regulation) REQUIRE r.regulation_id IS UNIQUE;",
        f"CREATE CONSTRAINT FOR (g:{cancer}_Gene) REQUIRE g.gene_id IS UNIQUE;",
    ]
    return query_list


def get_patient_query(cancer):
    return (f'LOAD CSV WITH HEADERS FROM "file:///Patient.csv" AS row\n'
            f"CREATE (:{cancer}_Patient:Patient {{patient_id: row.patient_id, mapping_id: row.mapping_id}});")


def get_gene_query(cancer):
    return (
            f'LOAD CSV WITH HEADERS FROM "file:///Gene.csv" AS row\n'
            'CALL {\n'
            'with row\n'
            f"CREATE (:{cancer}_Gene:Gene {{gene_id: row.gene_id, methylation: toFloat(row.methylation), mutation: toFloat(row.mutation)}})\n"
            "} IN TRANSACTIONS OF 5000 ROWS"
    )

def get_regulation_query(cancer):
    return (
            f'LOAD CSV WITH HEADERS FROM "file:///Regulation.csv" AS row\n'
            'CALL {\n'
            'with row\n'
            f"MATCH (source:{cancer}_Gene {{gene_id: row.source}})\n"
            f"MATCH (target:{cancer}_Gene {{gene_id: row.target}})\n"
            f"CREATE (source) -[:REGULATES]-> (:{cancer}_Regulation:Regulation {{regulation_id: row.regulation_id, fraction: toFloat(row.fraction), direction: row.direction, shared: toBoolean(row.shared)}}) -[:REGULATED]-> (target)\n"
            "} IN TRANSACTIONS OF 5000 ROWS")


def get_dysregulated_query(cancer):
    return (
            f'LOAD CSV WITH HEADERS FROM "file:///DYSREGULATED.csv" AS row\n'
            'CALL {\n'
            'with row\n'
            f"MATCH (patient:{cancer}_Patient {{mapping_id: row.mapping_id}})\n"
            f"MATCH (regulation:{cancer}_Regulation {{regulation_id: row.regulation_id}})\n"
            f"CREATE (patient) -[:DYSREGULATED {{value: row.value}}]-> (regulation)\n"
            "} IN TRANSACTIONS OF 5000 ROWS")

def get_methylated_query(cancer):
    return (
            f'LOAD CSV WITH HEADERS FROM "file:///METHYLATED.csv" AS row\n'
            'CALL {\n'
            'with row\n'
            f"MATCH (patient:{cancer}_Patient {{patient_id: row.patient_id}})\n"
            f"MATCH (gene:{cancer}_Gene {{gene_id: row.gene_id}})\n"
            f"CREATE (patient) -[:METHYLATED {{methylation: row.methylation}}]-> (gene)\n"
            "} IN TRANSACTIONS OF 5000 ROWS")

def get_mutated_query(cancer):
    return (
            f'LOAD CSV WITH HEADERS FROM "file:///MUTATED.csv" AS row\n'
            'CALL {\n'
            'with row\n'
            f"MATCH (patient:{cancer}_Patient {{patient_id: row.patient_id}})\n"
            f"MATCH (gene:{cancer}_Gene {{gene_id: row.gene_id}})\n"
            f"CREATE (patient) -[:MUTATED]-> (gene)\n"
            "} IN TRANSACTIONS OF 5000 ROWS")


def commit(query):
    db_connection = GraphDatabase.driver(uri=uri, auth=auth)
    session = db_connection.session()
    session.run(query)
    session.close()
    db_connection.close()


def commit_list(query_list):
    db_connection = GraphDatabase.driver(uri=uri, auth=auth)
    session = db_connection.session()
    for query in query_list:
        session.run(query)
    session.close()
    db_connection.close()

def write_csvs(network_path, stats_path):
    print("Load edge types")
    stats = pd.read_csv(stats_path, index_col=0)
    edge_types = dict(zip(stats.index.values, stats.coef_TF.values))


    print('Load feather')
    data = pd.read_feather(network_path)
    rownames = [str(patient_id) for patient_id in data.index.values]
    patient_ids = { key:value for key, value in zip(rownames, data["patient id"].tolist())}

    map = {}
    colnames = data.columns.values[1:]

    for colname in colnames:
        map[colname] = {key: value for key, value in zip(rownames, data[colname].tolist())}

    data = map


    print('Creating Patient.csv')
    write_patient_csv(patient_ids)

    print('Creating Regulation.csv')
    write_regulation_csv(data, edge_types)
    edge_types = None

    print('Creating DYSREGULATED.csv')
    write_dysregulated_csv(data)

    gene_set = get_gene_ids(data.keys())
    data = None

    print("Loading methylation")
    methylation_header = pd.read_csv(methylation_path, skipinitialspace=True, nrows=0).columns.tolist()
    columns = list(set(methylation_header) & set(patient_ids.values())) + ["Unnamed: 0"]
    methylation_data = pd.read_csv(methylation_path, skipinitialspace=True, usecols=columns)
    methylation_data.index = list(methylation_data.pop("Unnamed: 0"))
    methylation_means = dict(zip(methylation_data.index, list(methylation_data.mean(axis=1))))

    print("Creating METHYLATED.csv")
    write_methylation_csv(methylation_data)
    methylation_data = None

    print("Loading mutation")
    mutation_data = pd.read_csv(mutation_path, skipinitialspace=True,
                                usecols=["sample", "gene"])
    mutation_data = mutation_data[mutation_data['sample'].isin(set(patient_ids.values()))]
    mutation_patients_number = len(set(list(mutation_data["sample"])) & set(patient_ids.values()))
    mutation_pairs = list(set(zip(list(mutation_data["gene"]), list(mutation_data["sample"]))))
    mutation_data = None

    print("Creating MUTATED.csv")
    write_mutation_csv(mutation_pairs)

    gene_mutation_counts = Counter([gene for gene, patient in mutation_pairs])
    gene_mutation_frequency = {gene: count / mutation_patients_number for gene, count in gene_mutation_counts.items()}

    print('Creating Gene.csv')
    write_gene_csv(gene_set, methylation_means, gene_mutation_frequency)


def load_db(cancer):

    start = time.time()

    print("Committing cancer")
    commit(f"MERGE (:Cancer {{cancer_id: '{cancer}'}})")

    print("Committing pre-querys")
    commit_list(get_pre_query_list_cancer(cancer))
    print(time.time() - start)

    print("Committing patients")
    commit(get_patient_query(cancer))
    print(time.time() - start)

    print("Committing genes")
    commit(get_gene_query(cancer))
    print(time.time() - start)

    print("Committing regulations")
    commit(get_regulation_query(cancer))
    print(time.time() - start)

    print("Committing dysregulations")
    commit(get_dysregulated_query(cancer))
    print(time.time() - start)

    print("Committing methylation")
    commit(get_methylated_query(cancer))
    print(time.time() - start)

    print("Committing mutations")
    commit(get_mutated_query(cancer))
    print(time.time() - start)

    print("Deleting files")

    files = glob.glob(csv_path+'*')
    for f in files:
        os.remove(f)


if __name__ == '__main__':
    index = 0
    print("Committing overall pre-querys")
    commit_list(get_pre_query_list())

    # get cancer list from file names
    cancer_list = [network_path.split("/")[-3] for network_path in network_paths]


    for cancer, network_path, stats_path in zip(cancer_list, network_paths, stats_paths):
        print(cancer)
        write_csvs(network_path=network_path, stats_path=stats_path)
        load_db(cancer)


