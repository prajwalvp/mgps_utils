import argparse
import glob
import os


def get_tarballs(tarballs):
    complete_list = []
    for tarball in tarballs:
        complete_list.extend(glob.glob(tarball))
    return complete_list


def split_by_user(users, tarballs):
    nusers = len(users)
    ntarballs = len(tarballs)
    allocations = {i: [] for i in users}
    while tarballs:
        for user in users:
            if tarballs:
                allocations[user].append(tarballs.pop(0))
    return allocations


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tarballs", type=str, nargs="+", dest="tarballs", help="Tarballs to split [uses wildcards]")
    parser.add_argument("-u", "--users", type=str, nargs="+", dest="users", help="Users who will do the viewing")
    args = parser.parse_args()
    tarballs = get_tarballs(args.tarballs)
    allocations = split_by_user(args.users, tarballs)
    for user, allocation in allocations.items():
        print("{}: {}".format(user, " ".join(allocation)))

if __name__ == "__main__":
    main()