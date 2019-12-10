# import pyvxl, recursively importing all available submodules
def __import():

  import importlib
  import pkgutil

  # helper function: recursively import submodules
  def import_submodules(package):
    results = {}
    for loader, name, is_pkg in pkgutil.walk_packages(package.__path__):
      if not is_pkg: continue

      full_name = package.__name__ + '.' + name
      _module = importlib.import_module(full_name)
      results[full_name] = _module
      results.update(import_submodules(_module))

    return results

  # import all submodules
  main_module = importlib.import_module(__name__)
  modules = import_submodules(main_module)
  return modules


# save complete list of submodules
__sub_mod__ = list(__import().keys())

# list of top-level submodules
__all__ = [key.split('.')[-1] for key in __sub_mod__
           if key.count('.') == 1]
