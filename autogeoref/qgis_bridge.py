"""Run PyQGIS inside the user's open QGIS over the qgis_mcp socket. Never import qgis in-process."""
import json
import socket
import struct
import sys

HOST, PORT = "127.0.0.1", 9876
CLIENT_SITE = r"C:\Users\FAI-Akash\AppData\Local\uv\cache\archive-v0\6UbEkjU0r6KaKsVQ\Lib\site-packages"


class QgisUnavailable(RuntimeError):
    pass


def available(timeout=2.0):
    try:
        with socket.create_connection((HOST, PORT), timeout=timeout):
            return True
    except OSError:
        return False


def _send(sock, payload):
    data = json.dumps(payload).encode("utf-8")
    sock.sendall(struct.pack(">I", len(data)) + data)


def _recv(sock):
    head = b""
    while len(head) < 4:
        chunk = sock.recv(4 - len(head))
        if not chunk:
            raise QgisUnavailable("QGIS closed the connection")
        head += chunk
    n = struct.unpack(">I", head)[0]
    body = b""
    while len(body) < n:
        chunk = sock.recv(min(65536, n - len(body)))
        if not chunk:
            raise QgisUnavailable("QGIS closed the connection")
        body += chunk
    return json.loads(body.decode("utf-8"))


def execute(code, timeout=170):
    """Run `code` in QGIS. Returns {'ok': bool, 'stdout': str, 'error': str|None}."""
    if CLIENT_SITE not in sys.path:
        sys.path.append(CLIENT_SITE)
    try:
        from qgis_mcp.client import QgisMCPClient
        client = QgisMCPClient()
        client.connect()
        try:
            reply = client.send_command("execute_code", {"code": code, "timeout": timeout},
                                        timeout=timeout + 10)
        finally:
            client.disconnect()
    except ImportError:
        try:
            with socket.create_connection((HOST, PORT), timeout=timeout + 10) as sock:
                _send(sock, {"type": "execute_code", "params": {"code": code, "timeout": timeout}})
                reply = _recv(sock)
        except OSError as exc:
            raise QgisUnavailable(str(exc))
    result = reply.get("result", reply)
    return {"ok": bool(result.get("executed")), "stdout": result.get("stdout", ""),
            "error": result.get("error") or result.get("message")}


def json_result(code, timeout=170):
    """Run code whose last stdout line is JSON, and return the parsed object."""
    out = execute(code, timeout=timeout)
    if not out["ok"]:
        raise QgisUnavailable(out["error"] or "QGIS refused the code")
    lines = [ln for ln in out["stdout"].splitlines() if ln.strip()]
    if not lines:
        raise QgisUnavailable("QGIS produced no output")
    return json.loads(lines[-1])
