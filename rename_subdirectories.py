"""Rename the given subdirectory for each hall and AD."""

import argparse
import os

import common

def main(old_name_template, new_name_template):
    for site, ad in common.all_ads:
        os.rename(
            old_name_template.format(site=site, ad=ad),
            new_name_template.format(site=site, ad=ad),
        )
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('old_name_template')
    parser.add_argument('new_name_template')
    args = parser.parse_args()
    main(args.old_name_template, args.new_name_template)
