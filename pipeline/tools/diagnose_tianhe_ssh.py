#!/usr/bin/env python3
import os
import sys
import time

import pexpect


PASSWORD = os.environ.get("TIANHE_PASSWORD")
USER_HOST = "mdm@tianhe"
TIMEOUT = 20


def run_case(name, args):
    print(f"=== {name} ===", flush=True)
    print("argv:", " ".join(args), flush=True)
    child = pexpect.spawn(args[0], args[1:], encoding="utf-8", timeout=TIMEOUT)
    transcript = []
    sent_password = False
    start = time.time()
    status = "unknown"
    rc = None
    try:
        while True:
            idx = child.expect([
                r"(?i)password:",
                r"(?i)permission denied",
                r"(?i)are you sure you want to continue connecting",
                r"__OK__",
                pexpect.EOF,
                pexpect.TIMEOUT,
            ])
            transcript.append(child.before)
            if idx == 0:
                print("event=password_prompt", flush=True)
                if sent_password:
                    status = "second_password_prompt"
                    child.close(force=True)
                    break
                child.sendline(PASSWORD)
                sent_password = True
            elif idx == 1:
                status = "permission_denied"
                child.close(force=True)
                break
            elif idx == 2:
                print("event=hostkey_prompt", flush=True)
                child.sendline("yes")
            elif idx == 3:
                print("event=command_output_seen", flush=True)
                status = "saw_ok_waiting_eof"
            elif idx == 4:
                rc = child.exitstatus
                if status == "saw_ok_waiting_eof":
                    status = "success"
                else:
                    status = "eof_without_ok"
                break
            else:
                status = "timeout"
                child.close(force=True)
                break
    finally:
        if child.isalive():
            child.close(force=True)
        if rc is None:
            rc = child.exitstatus
    elapsed = time.time() - start
    text = "".join(transcript)
    print(f"status={status} rc={rc} elapsed={elapsed:.2f}s", flush=True)
    if text.strip():
        print("--- transcript tail ---", flush=True)
        print(text[-2000:], flush=True)
    print(flush=True)
    return status


def main():
    if not PASSWORD:
        print("TIANHE_PASSWORD is required for password-authentication diagnostics.", file=sys.stderr)
        return 2

    base_password = [
        "ssh",
        "-o", "GSSAPIAuthentication=no",
        "-o", "PreferredAuthentications=password",
        "-o", "PubkeyAuthentication=no",
        "-o", "NumberOfPasswordPrompts=1",
        "-o", "ConnectTimeout=10",
    ]
    cases = [
        ("password_no_tty_T", base_password + ["-T", USER_HOST, "echo __OK__"]),
        ("password_force_tty_tt", base_password + ["-tt", USER_HOST, "echo __OK__"]),
        ("password_bash_lc_no_tty", base_password + ["-T", USER_HOST, "bash", "-lc", "echo __OK__"]),
        ("default_batch_key_only", [
            "ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
            "-T", USER_HOST, "echo __OK__"
        ]),
    ]
    results = []
    for name, args in cases:
        results.append((name, run_case(name, args)))
    print("=== summary ===")
    for name, status in results:
        print(f"{name}\t{status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
