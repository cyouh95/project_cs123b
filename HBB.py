from utils import format_query_string, data_dir
import pandas as pd
import requests
import time
import json
import re
import os


clustalo_api = 'https://www.ebi.ac.uk/Tools/services/rest/clustalo'

clustalo_run_endpoint = clustalo_api + '/run'
clustalo_result_endpoint = clustalo_api + '/result'
clustalo_status_endpoint = clustalo_api + '/status'


class HBB:
    def __init__(self, basic_dataset, related_dataset):
        """
        Runs Clustal Omega analysis on input datasets and creates training and testing FASTA files for HMM
        """

        # Input is FASTA sequences downloaded from Blastp for basic and related organisms
        self.datasets = [basic_dataset, related_dataset]
        self.filenames = ['basic_dataset', 'related_dataset']
        self.seq_names = {}

        # Keep top 25 HBB sequences from each dataset to use as sample for HMM
        # Use 80% of that for training and 20% for testing
        self.gene = 'Hemoglobin subunit beta'
        self.num_seq = 25
        self.pct_training = 0.8

        # Keep 100bp segment from alignment to use for HMM
        # Define lowly conserved (ie. diverse) positions as having at least 5 different amino acids
        # Want segment to have between 10-15% lowly conserved positions
        self.window_size = 100
        self.diversity_min_aa = 5
        self.diversity_lower_bound = 0.1 * self.window_size
        self.diversity_upper_bound = 0.15 * self.window_size

        # Combine basic and related datasets to run MSA and phylogenetic tree analysis on Clustal Omega
        combined_cleaned_dataset = self.create_sample_dataset(self.seq_names)

        job_id = self.run_clustalo(combined_cleaned_dataset)

        while self.check_clustalo_job_status(job_id) != 'FINISHED':
            time.sleep(2)

        self.save_clustalo_phylotree_results(job_id)
        self.save_clustalo_alignment_results(job_id)

        # Determine best segment of alignment to keep
        self.alignment = self.parse_alignment()
        self.optimal_segment = self.get_optimal_alignment_segment()

        # Save training and testing datasets and alignment info
        self.gen_output_files()

        info = {
            'seq_names': self.seq_names,
            'optimal_segment': self.optimal_segment
        }

        json_object = json.dumps(info, indent=2)

        with open(f'{data_dir}/dataset_info.json', mode='w') as file:
            file.write(json_object)

    def create_sample_dataset(self, seq_names):
        """
        Saves cleaned sequence names for each dataset and returns combined FASTA
        """

        combined_cleaned_dataset = ''

        for filename, dataset in zip(self.filenames, self.datasets):
            cleaned_dataset = self.clean_fasta(dataset)

            # Save sequence names of sample
            seq_names[filename] = re.findall('>(\S+)', cleaned_dataset)

            combined_cleaned_dataset += cleaned_dataset

        with open(f'{data_dir}/clustalo_input.fasta', mode='w') as file:
            file.write(combined_cleaned_dataset)

        return combined_cleaned_dataset

    def clean_fasta(self, dataset):
        """
        Filters FASTA to include only top HBB sequences and renames sequence headers
        """

        cleaned_dataset = ''
        matches = re.findall('>[^>]+Full=' + self.gene + '([^ ;]+)?[^>]+\[([^\]]+)\](\n[^>]+)', dataset)

        for match in matches[:self.num_seq]:
            gene_suffix = match[0]
            organism = match[1].replace(' ', '_')
            sequence = match[2]

            cleaned_dataset += f'>{organism}_HBB{gene_suffix}{sequence}'

        return cleaned_dataset

    def get_optimal_alignment_segment(self):
        """
        Determine optimal segment in alignment to use for HMM

        Criteria:
            - Must meet diversity requirement (ie. have specific proportion of lowly conserved positions)
            - Maximize total number of lowly conserved and highly conserved positions
        """

        index = list(self.alignment.keys())  # sequence names
        data = [list(seq) for seq in self.alignment.values()]  # sequences

        # Create dataframe where each row is a sequence and each column is a base position in the alignment
        df = pd.DataFrame(data, index=index)

        # Add row containing number of unique amino acids for each position
        for pos in range(len(df.columns)):
            col = list(df.iloc[:50, pos])
            num_unique = len(set([aa for aa in col if aa != '-']))
            df.loc['num_unique', pos] = num_unique

        # Save as CSV to display as MSA results
        df.to_csv(f'{data_dir}/clustalo_alignment.csv', header=False)

        alignment_length = len(df.columns)

        optimal_segment = {
            'start': 0,
            'end': self.window_size,
            'num_lowly_conserved': None,
            'num_highly_conserved': None
        }

        # Want to maximize total number of lowly conserved and highly conserved positions
        max_desired_columns = 0

        for start in range(0, alignment_length - self.window_size + 1):
            end = start + self.window_size
            window = df.iloc[:, start:end]

            # Low conservation (e.g., 5+ AA)
            num_unique = window.loc['num_unique', :]
            num_lowly_conserved = len(num_unique[num_unique >= self.diversity_min_aa])

            # High conservation (marked `:` by Clustal Omega)
            markings = window.loc['markings', :]
            num_highly_conserved = len(markings[markings == ':'])

            # Total number of lowly conserved and highly conserved positions
            num_desired_columns = num_lowly_conserved + num_highly_conserved

            potential = ''
            if self.diversity_lower_bound <= num_lowly_conserved <= self.diversity_upper_bound \
                    and num_desired_columns >= max_desired_columns:
                max_desired_columns = num_desired_columns
                optimal_segment['start'] = start
                optimal_segment['end'] = end
                optimal_segment['num_lowly_conserved'] = num_lowly_conserved
                optimal_segment['num_highly_conserved'] = num_highly_conserved
                potential = '*'

            print(start, num_lowly_conserved, num_highly_conserved, potential, num_highly_conserved)

        print(optimal_segment)

        return optimal_segment

    def gen_output_files(self):
        """
        Trims aligned sequences to keep optimal segment and saves as training and testing dataset
        """

        start = self.optimal_segment['start']
        end = self.optimal_segment['end']

        cutoff = self.num_seq * self.pct_training

        for filename in self.filenames:
            training_file = open(f'{data_dir}/{filename}_training.fasta', mode='w')
            testing_file = open(f'{data_dir}/{filename}_testing.fasta', mode='w')

            for idx, seq_name in enumerate(self.seq_names[filename]):
                sequence = self.alignment[seq_name]
                trimmed_sequence = sequence[start:end]

                fasta_entry = f'>{seq_name}\n{trimmed_sequence}\n'

                # Split sequences into training and testing files
                if idx < cutoff:
                    training_file.write(fasta_entry)
                else:
                    testing_file.write(fasta_entry)

            training_file.close()
            testing_file.close()

    @staticmethod
    def parse_alignment():
        """
        Parse alignment results file and returns dictionary
        """

        alignment = {}

        # For determining position of the conservation markings
        start = None
        end = None

        with open(f'{data_dir}/clustalo_alignment.txt') as file:
            # Skip first heading line
            file.readline()

            for row in file:
                # Skip empty lines
                if not row.strip():
                    continue

                match = re.search('(\S+)\s+([-A-Z]+)', row)

                # MSA
                if match:
                    sequence_name = match.group(1)
                    sequence = match.group(2)
                    start = match.start(2)
                    end = match.end(2)

                # Conservation markings
                else:
                    sequence_name = 'markings'
                    sequence = row[start:end]

                if sequence_name not in alignment:
                    alignment[sequence_name] = ''
                alignment[sequence_name] += sequence

        return alignment

    @staticmethod
    def run_clustalo(dataset):
        """
        Runs Clustal Omega and returns job ID
        """

        params = {
            'email': os.environ['EMAIL_ADDRESS'],
            'guidetreeout': 'true',
            'dismatout': 'false',
            'dealign': 'false',
            'mbed': 'true',
            'mbediteration': 'true',
            'iterations': '0',
            'gtiterations': '-1',
            'hmmiterations': '-1',
            'outfmt': 'clustal',
            'order': 'aligned',
            'stype': 'protein',
            'sequence': dataset
        }
        data = format_query_string(params)

        headers = {
            'Content-Type': 'application/x-www-form-urlencoded',
            'Accept': 'text/plain'
        }

        r = requests.post(clustalo_run_endpoint, data=data, headers=headers)
        job_id = r.text

        return job_id

    @staticmethod
    def check_clustalo_job_status(job_id):
        """
        Checks Clustal Omega job status
        """

        r = requests.get(clustalo_status_endpoint + '/' + job_id)
        job_status = r.text

        return job_status

    @staticmethod
    def save_clustalo_alignment_results(job_id):
        """
        Saves multiple sequence alignment results from Clustal Omega job
        """

        r = requests.get(clustalo_result_endpoint + '/' + job_id + '/aln-clustal')
        alignment_results = r.text

        with open(f'{data_dir}/clustalo_alignment.txt', mode='w') as f:
            f.write(alignment_results)

    def save_clustalo_phylotree_results(self, job_id):
        """
        Saves phylogenetic tree results from Clustal Omega job
        """

        r = requests.get(clustalo_result_endpoint + '/' + job_id + '/phylotree')
        phylotree_results = r.text

        with open(f'{data_dir}/clustalo_phylotree.txt', mode='w') as f:
            f.write(phylotree_results)

        phylotree_json = {}

        phylotree_json = self.get_phylotree_json(phylotree_results, phylotree_json)

        json_object = json.dumps(phylotree_json, indent=2)

        with open(f'{data_dir}/clustalo_phylotree.json', mode='w') as file:
            file.write(json_object)

    def get_phylotree_json(self, phylotree_results, phylotree_json):
        """
        Creates json version of phylotree results for dendrogram
        """

        branch_regex = '\([^\(\)]+\)'
        leaf_regex = '([-|/\w]+)\n?:([\d.]+)'

        branches = re.finditer(branch_regex, phylotree_results)
        if not re.search(branch_regex, phylotree_results):
            root = list(phylotree_json.keys())[0]
            return {'name': root, 'children': phylotree_json[root]}

        for branch in branches:
            leaves = re.finditer(leaf_regex, branch.group())
            leaf_names = []
            children = []

            for leaf in leaves:
                leaf_name = leaf.group(1)
                leaf_score = leaf.group(2)

                leaf_names.append(leaf_name)

                json_entry = {'name': leaf_name, 'score': leaf_score}
                children.append(json_entry)

                if leaf_name in phylotree_json:
                    json_entry['children'] = phylotree_json[leaf_name]
                    phylotree_json.pop(leaf_name)

            branch_name = '|'.join(leaf_names)

            phylotree_json[branch_name] = children
            phylotree_results = phylotree_results.replace(branch.group(), branch_name)

        return self.get_phylotree_json(phylotree_results, phylotree_json)
