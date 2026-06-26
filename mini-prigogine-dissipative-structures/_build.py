import os
BASE = r"F:\nano-everything\mini-complex-control-theory\26. mini-synergetics-dissipative\mini-prigogine-dissipative-structures"
def wf(relpath, content):
    path = os.path.join(BASE, relpath)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    return content.count("\n")
tot = 0