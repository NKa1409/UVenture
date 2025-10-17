import os
import shlex
import time
import paramiko

class NetworkTaskScheduler:
    """A class to manage network tasks."""
    
    def __init__(self):
        """Initialize the task scheduler with an empty task list."""
        self.tasks = []
        self.ssh_connections = []
        self.remote_path = os.path.normpath(os.path.join("./", 'remote_projects', "UVenture")) # Default remote path
        # Get the upper folder of this script on this computer
        self.local_path = os.path.normpath(os.path.join(os.getcwd(), '..'))  # Default local path

    def add_ssh_connection(self, host, port, username, password, key_file=None, infos=None):
        """Add a new SSH connection to the list."""
        connection = {
            'host': host,
            'port': port,
            'username': username,
            'password': password,
            'key_file': key_file,
            'infos': infos
        }
        self.ssh_connections.append(connection)

    def remove_ssh_connection(self, host, port):
        """Remove an SSH connection from the list."""
        self.ssh_connections = [conn for conn in self.ssh_connections if not (conn['host'] == host and conn['port'] == port)]

    def create_connection(self, host, port, username, password, key_file=None):
        """Create a new SSH connection."""
        client = paramiko.SSHClient()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        if key_file:
            client.connect(hostname=host, port=port, username=username, key_filename=key_file)
        else:
            client.connect(hostname=host, port=port, username=username, password=password)
        return client
    
    def close_ssh_connection(self, ssh_client):
        """Close an SSH connection."""
        ssh_client.close()
        return True

    def upload_project_folder(self, ssh_client, local_path):
        """Upload a project folder to the remote server."""
        sftp = ssh_client.open_sftp()
        sftp.put(local_path, self.remote_path)
        sftp.close()

    def delete_project_folder(self, ssh_client):
        """Delete a project folder on the remote server."""
        sftp = ssh_client.open_sftp()
        try:
            sftp.remove(self.remote_path)
        except FileNotFoundError:
            print(f"Remote path {self.remote_path} not found.")
        sftp.close()

    def start_program_on_remote(self, ssh_client, program_name):
        """Start a program on the remote server."""
        operating_system_type = self.get_operating_system_type(ssh_client)
        command = ''
        if operating_system_type == 'Linux':
            command = f'cd {self.remote_path} && source venv/bin/activate && nohup python3 {program_name} > output.log 2>&1 &'
        elif operating_system_type == 'Darwin':
            command = f'cd {self.remote_path} && source venv/bin/activate && nohup python3 {program_name} > output.log 2>&1 &'
        elif operating_system_type == 'Windows':
            command = f'START /B "{self.remote_path}\\venv\\Scripts\\python.exe" "{self.remote_path}\\{program_name}"'
        else:
            print(f"Unsupported operating system type: {operating_system_type}")
            return False
        bash_command = f"bash -c {shlex.quote(command)}" if operating_system_type in ['Linux', 'Darwin'] else command
        stdin, stdout, stderr = ssh_client.exec_command(bash_command)
        time.sleep(1)  # Wait for the command to execute
        stderr_output = stderr.read().decode()
        exit_status = stdout.channel.recv_exit_status()
        if exit_status == 0:
            print(f"Program {program_name} started successfully on remote server.")
            return True
        else:
            print(f"Failed to start program {program_name} on remote server. Error: {stderr_output}")
            return False

    def get_core_count_on_remote(self, ssh_client):
        """Get the number of CPU cores on the remote server."""
        operating_system_type = self.get_operating_system_type(ssh_client)
        command = ''
        if operating_system_type == 'Linux':
            command = 'nproc'
        elif operating_system_type == 'Darwin':
            command = 'sysctl -n hw.ncpu'
        elif operating_system_type == 'Windows':
            command = 'wmic cpu get NumberOfCores'
        else:
            print(f"Unsupported operating system type: {operating_system_type}")
            return None
        
        stdin, stdout, stderr = ssh_client.exec_command(command)
        core_count = stdout.read().decode().strip()
        if core_count.isdigit():
            return int(core_count)
        else:
            print(f"Error retrieving core count: {stderr.read().decode()}")
            return None

    def get_operating_system_type(self, ssh_client):
        """Get the operating system type of the remote server."""
        stdin, stdout, stderr = ssh_client.exec_command('uname -s')
        os_type = stdout.read().decode().strip()
        if os_type in ['Linux', 'Darwin']:
            return os_type
        else:
            if stderr.read().decode():
                print(f"Error determining OS type: {stderr.read().decode()}")
                print("Defaulting to Windows-like behavior.")
            return 'Windows'






if __name__ == "__main__":
    network_scheduler = NetworkTaskScheduler()
    print(network_scheduler.local_path)
    print(network_scheduler.remote_path)
