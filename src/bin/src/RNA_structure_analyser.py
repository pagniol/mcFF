import subprocess
import paramiko
import os

from utils.extract_bps import *

class RNASecondaryStructureAnalyzer:
    def __init__(self, rnaview_path="rnaview", x3dna_path="x3dna-dssr",
                 mcannotate_remote_path="./stage-E24/tools/MC-Annotate",
                 ssh_host="sagniol@cluster.iric.ca"):
        """
        Initialise l'analyseur pour générer des structures dot-bracket.

        :param rnaview_path: Chemin vers l'exécutable RNAView (local)
        :param x3dna_path: Chemin vers l'exécutable x3DNA DSSR (local)
        :param mcannotate_remote_path: Chemin distant vers MC-Annotate
        :param ssh_host: Hôte SSH pour exécuter MC-Annotate
        """
        self.rnaview_path = rnaview_path
        self.x3dna_path = x3dna_path
        self.mcannotate_remote_path = mcannotate_remote_path
        self.ssh_host = ssh_host
        
    def run_rnaview(self, pdb_file):
        """
        Exécute RNAView localement sur un fichier PDB et retourne la sortie.

        :param pdb_file: Chemin vers le fichier PDB
        :return: Sortie RNAView sous forme de chaîne
        """
        try:
            result = subprocess.run([self.rnaview_path, pdb_file],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            return result.stdout.decode('utf-8')
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Erreur lors de l'exécution de RNAView : {e.stderr.decode('utf-8')}")

    def run_x3dna(self, pdb_file):
        """
        Exécute x3DNA DSSR localement sur un fichier PDB et retourne la sortie.

        :param pdb_file: Chemin vers le fichier PDB
        :return: Sortie x3DNA DSSR sous forme de chaîne
        """
        try:
            result = subprocess.run([self.x3dna_path, "--input", pdb_file],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            return result.stdout.decode('utf-8')
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Erreur lors de l'exécution de x3DNA DSSR : {e.stderr.decode('utf-8')}")

    def run_mcannotate(self, pdb_file):
        """
        Exécute MC-Annotate sur un serveur distant via SSH et retourne la sortie.

        :param pdb_file: Chemin vers le fichier PDB
        :return: Sortie MC-Annotate sous forme de chaîne
        """
        remote_pdb_file = f"/tmp/{os.path.basename(pdb_file)}"
        ssh_client = paramiko.SSHClient()
        ssh_client.load_system_host_keys()
        ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        try:
            # Connexion au serveur
            ssh_client.connect(self.ssh_host.split('@')[1], username=self.ssh_host.split('@')[0])

            # Copier le fichier PDB sur le serveur
            with paramiko.SFTPClient.from_transport(ssh_client.get_transport()) as sftp:
                sftp.put(pdb_file, remote_pdb_file)

            # Exécuter MC-Annotate
            command = f"{self.mcannotate_remote_path} {remote_pdb_file}"
            stdin, stdout, stderr = ssh_client.exec_command(command)

            # Collecter les sorties
            output = stdout.read().decode('utf-8')
            error_message = stderr.read().decode('utf-8')
            if error_message:
                raise RuntimeError(f"Erreur MC-Annotate : {error_message}")

            # Supprimer le fichier distant
            ssh_client.exec_command(f"rm {remote_pdb_file}")

        finally:
            ssh_client.close()

        return output

    def generate_dot_bracket(self, pdb_file):
        """
        Génère une structure dot-bracket en combinant les résultats de RNAView, x3DNA DSSR et MC-Annotate.

        :param pdb_file: Chemin vers le fichier PDB
        :return: Dot-bracket
        """
        # Exécuter les outils et récupérer les sorties
        rnaview_output = self.run_rnaview(pdb_file)
        x3dna_output = self.run_x3dna(pdb_file)
        mcannotate_output = self.run_mcannotate(pdb_file)

        # Extraire les paires de bases de chaque outil
        rnaview_pairs = extract_rnaview_bps(rnaview_output.splitlines())
        x3dna_pairs = extract_x3dna_bps(x3dna_output.splitlines())
        mcannotate_pairs = extract_mcannotate_bps(mcannotate_output.splitlines())

        # Combiner les paires et générer le dot-bracket
        combined_pairs = combine_pairs(rnaview_pairs + x3dna_pairs + mcannotate_pairs)
        sequence_length = len(extract_sequence_from_pdb(pdb_file))
        return generate_dot_bracket(combined_pairs, sequence_length)

# Exemple d'utilisation
if __name__ == "__main__":
    analyzer = RNASecondaryStructureAnalyzer()
    pdb_file = "example.pdb"
    dot_bracket = analyzer.generate_dot_bracket(pdb_file)
    print(f"Dot-bracket:\n{dot_bracket}")
