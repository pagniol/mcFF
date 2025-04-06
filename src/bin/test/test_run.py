import subprocess

def execute_command(command, capture_output=True, shell=False):
    """
    Exécute une commande système et retourne le résultat.

    :param command: Liste des arguments de la commande ou chaîne de commande complète si shell=True
    :param capture_output: Si True, capture stdout et stderr
    :param shell: Si True, exécute la commande via le shell
    :return: Tuple (stdout, stderr, returncode)
    :raises: subprocess.CalledProcessError si la commande échoue
    """
    try:
        result = subprocess.run(
            command,
            capture_output=capture_output,
            text=True,
            shell=shell,
            check=True
        )
        return result.stdout, result.stderr, result.returncode
    except subprocess.CalledProcessError as e:
        # Retourne les informations en cas d'échec
        return e.stdout, e.stderr, e.returncode
    
output, error, code = execute_command(["/u/sagnioln/stage-E24/tools/MC-Annotate", "/u/sagnioln/stage-E24/test/6ugg.pdb"])
print("Sortie standard :", output)
print("Sortie d'erreur :", error)
print("Code de retour :", code)
