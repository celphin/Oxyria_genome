from pathlib import Path
directory = Path('.').absolute()
folders = directory.glob('dir.*')

for f in folders:
    split = f.as_posix().split(".")
    dir_name = "dir" "." + str(int(split[1])).zfill(1)
    f.rename(f.with_name(dir_name))
    print(dir_name)
