def dict_contains(check_key, dict):
    for key in dict.keys():
        if key == check_key:
            return True
    return False