import os
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


def copy_tarballs(user, localpath, basepath, tarballs):
    ssh = paramiko.SSHClient()
    ssh.load_host_keys(os.path.expanduser(
        os.path.join("~", ".ssh", "known_hosts")))
    ssh.connect(SERVER,
                username=user,
                password=getpass.getpass(prompt="Password: "))
    sftp = ssh.open_sftp()
    for tarball in tarballs:
        print(os.path.join(basepath, tarball),
              os.path.join(localpath, tarball))
        sftp.get(os.path.join(basepath, tarball),
                 os.path.join(localpath, tarball),
                 callback=progress)
    print("\nCopy complete")
    sftp.close()
    ssh.close()


def is_nonempty_tar_file(archive):
    with tarfile.open(archive, "r") as tar:
        try:
            file_content = tar.getmembers()
            return len(file_content) > 0
        except Exception as exc:
            print(f"Reading tarfile failed for {archive}")


def unpack_merge(localpath, tarballs):
    merged_cands_filename = os.path.join(localpath, "merged_candidates.csv")
    merged_cands_file = open(merged_cands_filename, "wb")
    cands_file = os.path.join(localpath, "candidates.csv")
    has_header = False

    for tarball in tarballs:
        print(tarball)
        f = tarfile.open(os.path.join(localpath, tarball), "r:gz")
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
    shutil.move(merged_cands_filename, cands_file)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--user", type=str, dest="user",
                        help="Username for SSH connection to MPIfR", default=getpass.getuser())
    parser.add_argument("-r", "--remotepath", type=str, dest="remotepath",
                        help="Remote path where tarballs are located")
    parser.add_argument("-l", "--localpath", type=str, dest="localpath",
                        help="Localpath where tarballs should be unpacked and merged", default="./")
    parser.add_argument("-t", "--tarballs", type=str, dest="tarballs",
                        nargs="+", help="List of tarballs to retrieve and extract")
    parser.add_argument("-s", "--skip_dload",  type=int, dest="skip_dload_flag",
                        help="Skip download if tarballs are already downloaded (Default=0)", default=0)
    args = parser.parse_args()

    if not args.skip_dload_flag:
        copy_tarballs(args.user, args.localpath,
                      args.remotepath, args.tarballs)
    unpack_merge(args.localpath, args.tarballs)


if __name__ == "__main__":
    main()
