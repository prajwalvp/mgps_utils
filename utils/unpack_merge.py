import os
import sys
import paramiko
import getpass
import tarfile
import shutil
import fcntl
import termios
import struct
import argparse

SERVER = "portal.mpifr-bonn.mpg.de"

def terminal_size():
    h, w, hp, wp = struct.unpack('HHHH',
        fcntl.ioctl(0, termios.TIOCGWINSZ,
        struct.pack('HHHH', 0, 0, 0, 0)))
    return w, h

def progress(transferred, total):
    w, _ = terminal_size()
    w = w//2
    complete = int(w * transferred/total)
    print("Copying: |{}>{}| {}%".format("="*complete, " "*(w-complete),
        int(100*(transferred/total))), end="\r")


def connect(user):
    ssh = paramiko.SSHClient()
    ssh.load_host_keys(os.path.expanduser(os.path.join("~", ".ssh", "known_hosts")))
    ssh.connect("portal.mpifr-bonn.mpg.de",
        username=user,
        password=getpass.getpass(prompt="Password: "))
    return ssh


def find_tarballs(conn, basepath, ids):
    print("Finding matching tarballs")
    stdin, stdout, stderr = conn.exec_command(f"ls -1 {basepath}")
    tarballs = []
    for fname in stdout.readlines():
        fname = fname.strip()
        if int(fname.split("_")[1]) in ids:
            print(f"Found: {fname}")
            tarballs.append(fname)
    return tarballs


def copy_tarballs(conn, localpath, basepath, tarballs):
    sftp = conn.open_sftp()
    for tarball in tarballs:
        print(os.path.join(basepath, tarball), os.path.join(localpath, tarball))
        sftp.get(os.path.join(basepath, tarball),
            os.path.join(localpath, tarball),
            callback=progress)
    print("\nCopy complete")
    sftp.close()


def unpack_merge(localpath, tarballs):
    merged_cands_file = open(os.path.join(localpath, "merged_candidates.csv"), "wb")
    has_header = False

    for tarball in tarballs:
        f = tarfile.open(localpath+'/'+tarball, "r:gz")
        candfile = f.extractfile("candidates.csv")
        if not has_header:
            merged_cands_file.write(candfile.read())
            has_header = True
        else:
            candfile.readline()
            for line in candfile.readlines():
                merged_cands_file.write(line)
        f.extractall(path=localpath)
        f.close()
    merged_cands_file.close()
    shutil.move("{}/merged_candidates.csv".format(localpath), "{}/candidates.csv".format(localpath))


# Remember to close the conn

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--user", type=str, dest="user", help="Username for SSH connection to MPIfR", default=getpass.getuser())
    parser.add_argument("-r", "--remotepath", type=str, dest="remotepath", help="Remote path where tarballs are located")
    parser.add_argument("-l", "--localpath", type=str, dest="localpath", help="Localpath where tarballs should be unpacked and merged", default="./")
    parser.add_argument("-s", "--start", type=int, dest="start", help="First pointing ID to retrieve")
    parser.add_argument("-e", "--end", type=int, dest="end", help="Last pointing ID to retrieve")
    args = parser.parse_args()
    idxs = list(range(args.start, args.end+1))
    conn = connect(args.user)
    tarballs = find_tarballs(conn, args.remotepath, idxs)
    copy_tarballs(conn, args.localpath, args.remotepath, tarballs)
    conn.close()
    unpack_merge(args.localpath, tarballs)

if __name__ == "__main__":
    main()

