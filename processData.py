from sys import argv
from os import path
from struct import pack

if __name__ == "__main__":
  in_dir = argv[1]
  out_dir = argv[2]
  param_file = open(path.join(out_dir, "params"), "w+b")
  with open(path.join(in_dir, "ep.txt"), "r") as inf:
    param_file.write(int(inf.read()).to_bytes(4, 'little'))
  with open(path.join(in_dir, "max_level.txt"), "r") as inf:
    M = int(inf.read())
    param_file.write(M.to_bytes(4, 'little'))

  with open(path.join(out_dir, "index"), "w+b") as index:
    with open(path.join(in_dir, "index.txt"), "r") as inf:
      index.writelines(list(map(lambda x: int(x).to_bytes(4, 'little', signed=True), inf.readlines())))

  with open(path.join(out_dir, "indptr"), "w+b") as indptr:
    with open(path.join(in_dir, "indptr.txt"), "r") as inf:
      lines = inf.readlines()
      l = int(len(lines)) - 1
      param_file.write(l.to_bytes(4, 'little'))
      indptr.writelines(list(map(lambda x: int(x).to_bytes(4, 'little'), lines)))

  with open(path.join(out_dir, "level_offset"), "w+b") as level_offset:
    with open(path.join(in_dir, "level_offset.txt"), "r") as inf:
      level_offset.writelines(list(map(lambda x: int(x).to_bytes(4, 'little'), inf.readlines())))

  with open(path.join(out_dir, "vect"), "w+b") as vect:
    with open(path.join(in_dir, "vect.txt"), "r") as inf:
      lines = inf.readlines()
      d = int(len(lines[0].split()))
      param_file.write(d.to_bytes(4, 'little'))
      for line in lines:
        vect.writelines(list(map(lambda x: pack("<f", float(x)), line.split())))
