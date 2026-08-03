import os


def hmm_profiles_in_dirs(directories):
    profiles = []
    for directory in directories or []:
        directory = os.path.abspath(directory)
        if not os.path.isdir(directory):
            raise ValueError(f"HMM profile directory does not exist: {directory}")
        profiles.extend(
            os.path.join(directory, name)
            for name in sorted(os.listdir(directory))
            if name.endswith(".hmm") and os.path.isfile(os.path.join(directory, name))
        )
    return profiles
