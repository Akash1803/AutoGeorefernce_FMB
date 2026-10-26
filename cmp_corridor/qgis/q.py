import json, socket, struct, sys
def call(cmd, params=None, timeout=120):
    s = socket.create_connection(("127.0.0.1", 9876), timeout=timeout)
    msg = json.dumps({"type": cmd, "params": params or {}}).encode()
    s.sendall(struct.pack(">I", len(msg)) + msg)
    hdr = b""
    while len(hdr) < 4: hdr += s.recv(4 - len(hdr))
    n = struct.unpack(">I", hdr)[0]; buf = b""
    while len(buf) < n: buf += s.recv(min(65536, n - len(buf)))
    s.close(); return json.loads(buf)
def code(src, timeout=300):
    r = call("execute_code", {"code": src}, timeout); return r
if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    src = sys.stdin.read()
    r = code(src)
    out = r.get("result", r)
    print(out if isinstance(out, str) else json.dumps(out, ensure_ascii=False, indent=1)[:6000])
