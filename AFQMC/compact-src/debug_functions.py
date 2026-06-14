import sys


def whos():
    # Get all variables in the local scope
    vars_dict = locals()
    
    # Calculate size and type for each variable
    print(f"{'Variable':<20} {'Type':<20} {'Size (Bytes)':<15}")
    print("-" * 55)
    
    for name, value in vars_dict.items():
        if not name.startswith('_'): # Filter out internal Python variables
            size = sys.getsizeof(value)
            type_name = type(value).__name__
            print(f"{name:<20} {type_name:<20} {size:<15}")