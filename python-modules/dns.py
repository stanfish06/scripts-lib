#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Check whether the current network's DNS redirects a domain to a block page.

Compares the system resolver's answer with a public DNS-over-HTTPS answer, then
connects to each address and prints the certificate it serves.

Usage: uv run check_dns_block.py [host]
Exit code: 0 = reachable with a valid certificate, 1 = blocked or failed.
"""

import _ssl
import argparse
import json
import socket
import ssl
import sys
import tempfile
import urllib.request
from pathlib import Path

# Domain to test (override with the first CLI argument)
HOST = "docs.scloop.dev"
PORT = 443
# Public DNS-over-HTTPS resolvers (JSON API), tried in order; 1.1.1.1 by IP so it skips system DNS
DOH_URLS = ["https://1.1.1.1/dns-query", "https://dns.google/resolve"]
# Seconds per network call
TIMEOUT = 10


def system_dns(host: str) -> list[str]:
    """IPv4 addresses from the OS resolver (the one the browser uses)."""
    infos = socket.getaddrinfo(host, PORT, socket.AF_INET, socket.SOCK_STREAM)
    return sorted({str(info[4][0]) for info in infos})


def public_dns(host: str) -> tuple[str, list[str]]:
    """IPv4 addresses from the first public DoH resolver that answers."""
    for url in DOH_URLS:
        req = urllib.request.Request(
            f"{url}?name={host}&type=A", headers={"accept": "application/dns-json"}
        )
        try:
            with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
                answer = json.load(resp).get("Answer", [])
        except (OSError, ValueError) as err:
            print(f"  {url} failed: {err}")
            continue
        # type 1 = A record; skip CNAME hops
        return url, sorted(a["data"] for a in answer if a["type"] == 1)
    return "", []


def served_cert(host: str, ip: str) -> dict:
    """Certificate the server at ip presents for host, read without verification."""
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    with (
        socket.create_connection((ip, PORT), timeout=TIMEOUT) as sock,
        ctx.wrap_socket(sock, server_hostname=host) as tls,
    ):
        der = tls.getpeercert(binary_form=True)
    if der is None:
        raise ssl.SSLError(f"{ip} sent no certificate")
    # getpeercert() returns {} without verification; decode the raw cert with CPython's internal helper
    with tempfile.TemporaryDirectory() as tmp:
        pem = Path(tmp) / "cert.pem"
        pem.write_text(ssl.DER_cert_to_PEM_cert(der))
        return _ssl._test_decode_cert(str(pem))  # pyright: ignore[reportAttributeAccessIssue]


def verified_head(host: str, ip: str) -> tuple[bool, str]:
    """HEAD / on ip with the browser's checks (trusted chain + hostname match)."""
    ctx = ssl.create_default_context()
    try:
        with (
            socket.create_connection((ip, PORT), timeout=TIMEOUT) as sock,
            ctx.wrap_socket(sock, server_hostname=host) as tls,
        ):
            tls.sendall(
                f"HEAD / HTTP/1.1\r\nHost: {host}\r\nConnection: close\r\n\r\n".encode()
            )
            return True, tls.recv(256).split(b"\r\n", 1)[0].decode()
    # Hostname mismatch or untrusted chain: what the browser reports as ERR_CERT_*
    except ssl.SSLCertVerificationError as err:
        return False, err.verify_message
    except OSError as err:
        return False, str(err)


def dn(rdns: tuple) -> str:
    return ", ".join(f"{k}={v}" for rdn in rdns for k, v in rdn)


def report(host: str, ip: str) -> bool:
    """Print the served certificate and the verified result for one address."""
    print(f"  {ip}")
    try:
        cert = served_cert(host, ip)
    except OSError as err:
        print(f"    connect failed: {err}")
        return False
    print(f"    subject: {dn(cert['subject'])}")
    print(f"    issuer:  {dn(cert['issuer'])}")
    print(f"    SANs:    {', '.join(v for _, v in cert.get('subjectAltName', ()))}")
    ok, detail = verified_head(host, ip)
    print(f"    verified TLS: {'OK' if ok else 'FAIL'}, {detail}")
    return ok


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("host", nargs="?", default=HOST)
    host = parser.parse_args().host

    print("== System DNS")
    try:
        sys_ips = system_dns(host)
    except socket.gaierror as err:
        print(f"  lookup failed: {err}")
        sys_ips = []
    print(f"  {', '.join(sys_ips) or 'none'}")

    print("== Public DNS")
    url, pub_ips = public_dns(host)
    print(f"  {url or 'all resolvers failed'}: {', '.join(pub_ips) or 'none'}")

    print("== Certificates via system DNS")
    sys_ok = [report(host, ip) for ip in sys_ips]

    pub_ok = False
    if pub_ips and set(pub_ips) != set(sys_ips):
        print("== Certificate via public DNS")
        pub_ok = report(host, pub_ips[0])

    print("== Result")
    if sys_ips and all(sys_ok):
        print(f"  OK: this network reaches {host} with a valid certificate.")
        return 0
    if pub_ips and not set(sys_ips) & set(pub_ips):
        print(
            f"  BLOCKED: this network's DNS answers {', '.join(sys_ips) or 'nothing'}"
            f" instead of {', '.join(pub_ips)}."
        )
        if pub_ok:
            print(
                "  The site itself is fine: the public-DNS address serves a valid certificate."
            )
        return 1
    print("  FAIL: TLS check failed; see the details above.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
