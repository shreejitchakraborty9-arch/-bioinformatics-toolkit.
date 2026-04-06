import os
with open(r'c:\essential files\python\01-setup\ANTIGRAVITY\bioinformatics-toolkit\app.py', 'r', encoding='utf-8') as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if 'deferred: {e}' in line:
        lines[i] = '        app.logger.warning(f"Database initialization deferred: {e}")\n'
        print(f"Fixed line {i+1}")

with open(r'c:\essential files\python\01-setup\ANTIGRAVITY\bioinformatics-toolkit\app.py', 'w', encoding='utf-8') as f:
    f.writelines(lines)
