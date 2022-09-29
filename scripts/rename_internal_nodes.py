import sys


def main(path_to_dist_tree, path_to_named_tree, path_to_out):
    # Iteratively going from leaves to root and add names to dist tree
    pass


if __name__ == "__main__":
    try:
        path_to_dist_tree, path_to_named_tree, path_to_out = sys.argv[1:]
        main(path_to_dist_tree, path_to_named_tree, path_to_out)
    except:
        print("ERROR\nUSAGE: script.py path_to_dist_tree path_to_named_tree path_to_out_tree", file=sys.stderr)
