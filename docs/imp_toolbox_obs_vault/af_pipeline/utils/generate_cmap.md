```python
def generate_cmap(n, scheme="soft-warm"):
    """ Generate a list of n colors
    (modfied from chatgpt)

    Args:
        n (int): number of colors
        scheme (str, optional): Defaults to "soft-warm".

    Returns:
        colors (list): list of n colors
    """

    import random
    import time

    colors = set()
    start = time.time()

    assert n > 0; "Number of colors must be greater than 0"
    if n != 2:
        assert scheme != "binary"; "Binary scheme is only valid for 2 colors"

    if n == 2 and scheme == "binary":
        return ["black", "green"]

    while len(colors) < n:
        if scheme == "standard":
            r = random.randint(0, 255)
            g = random.randint(0, 255)
            b = random.randint(0, 255)

        elif scheme == "non-bright":
            r = random.randint(50, 180)
            g = random.randint(50, 180)
            b = random.randint(50, 180)

        elif scheme == "earth-tone":
            r = random.randint(100, 180)
            g = random.randint(60, 140)
            b = random.randint(40, 120)

        elif scheme == "cool-tone":
            r = random.randint(50, 120)
            g = random.randint(100, 180)
            b = random.randint(120, 255)

        elif scheme == "soft-warm":
            r = random.randint(180, 255)
            g = random.randint(130, 200)
            b = random.randint(90, 160)

        elif scheme == "contrasting-non-bright":
            if len(colors) % 2 == 0:  # Alternate between darker and lighter muted tones
                r = random.randint(40, 120)
                g = random.randint(40, 120)
                b = random.randint(40, 120)
            else:
                r = random.randint(140, 200)
                g = random.randint(140, 200)
                b = random.randint(140, 200)

        else:
            raise ValueError(
                "Invalid scheme. Choose from 'non-bright', 'earth-tone', 'cool-tone', or 'soft-warm'."
            )

        color = "#{:02x}{:02x}{:02x}".format(r, g, b)
        colors.add(color)

        if time.time() - start > 10:
            break

    return list(colors)
```