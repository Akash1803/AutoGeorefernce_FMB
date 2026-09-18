import sys, json
sys.path.insert(0, r"C:\Users\FAI-Akash\AppData\Local\uv\cache\archive-v0\6UbEkjU0r6KaKsVQ\Lib\site-packages")
from qgis_mcp.client import QgisMCPClient
code = open(sys.argv[1], encoding="utf-8").read()
c = QgisMCPClient(); c.connect()
r = c.send_command("execute_code", {"code": code, "timeout": 170}, timeout=180)
c.disconnect()
res = r.get("result", r)
out = res.get("stdout", "")
err = res.get("stderr", "")
print("STATUS:", r.get("status"), "| executed:", res.get("executed"), "| error:", res.get("error") or res.get("message"))
if err: print("STDERR:", err[:2000])
try:
    print(json.dumps(json.loads(out.strip().splitlines()[-1]), indent=1))
except Exception:
    print("STDOUT:", out[:4000])
