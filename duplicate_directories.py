"""Create a new directory structure with many folders symlinked.

Specifically, link to existing versions of:
    - processed_ad?
    - hadded_ad?
    - singles_ad?
    - acc_ad?

And create the following empty directories:
    - sub_ad?

The directories are linked and created for both "nominal" and "adtime"
variants.
"""
import argparse
import os

import common


empty_dirs = [
    'EH{site}/sub_nominal_ad{ad}/',
    'EH{site}/sub_using_adtime_ad{ad}/',
]


linked_dirs = [
    'EH{site}/processed_ad{ad}',
    'EH{site}/hadded_ad{ad}',
    'EH{site}/singles_ad{ad}',
    'EH{site}/acc_nominal_ad{ad}',
    'EH{site}/acc_using_adtime_ad{ad}',
]

def main(base_location, output_base_location):
    os.makedirs(output_base_location, exist_ok=True)
    for site, ad in common.all_ads:
        for dirname_template in empty_dirs:
            pathname = os.path.join(
                output_base_location,
                dirname_template.format(site=site, ad=ad),
            )
            os.makedirs(pathname, exist_ok=True)
        for dirname_template in linked_dirs:
            dirname = dirname_template.format(site=site, ad=ad)
            os.symlink(
                os.path.join(base_location, dirname),
                os.path.join(output_base_location, dirname),
            )
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('base_location')
    parser.add_argument('output_base_location')
    args = parser.parse_args()
    main(args.base_location, args.output_base_location)
