How to import the yaml config file in a python script:

    import yaml
    path_i = os.path.join(os.environ["PROJ_irox"], "config", "config.yml")
    with open(path_i) as file:
        config_dict = yaml.load(file, Loader=yaml.FullLoader)

    api_key = config_dict["mpcontrib"]["api_key"]
