#!/usr/bin/env python3

from __future__ import annotations

import http.client
import ssl
import sys
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlsplit


LISTEN_HOST = "127.0.0.1"
LISTEN_PORT = 8765

UPSTREAM_HOST = "downloads.pacbcloud.com"
UPSTREAM_PREFIX = (
    "/public/dataset/"
    "HG002-CpG-methylation-202202/"
)

ALLOWED_FILES = {
    "HG002.GRCh38.haplotagged.bam",
    "HG002.GRCh38.haplotagged.bam.bai",
}

BUFFER_SIZE = 1024 * 1024
MAX_OPEN_ATTEMPTS = 12

TLS_CONTEXT = ssl.create_default_context()
TLS_CONTEXT.set_alpn_protocols(["http/1.1"])


def log(message: str) -> None:
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(
        f"[{timestamp}] {message}",
        file=sys.stderr,
        flush=True,
    )


class ProxyHandler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"
    server_version = "HG002HTTP1Proxy/1.0"

    def log_message(
        self,
        format_string: str,
        *args,
    ) -> None:
        log(
            f"{self.client_address[0]} "
            f"{format_string % args}"
        )

    def selected_file(self) -> str | None:
        filename = (
            urlsplit(self.path)
            .path
            .rsplit("/", 1)[-1]
        )

        if filename in ALLOWED_FILES:
            return filename

        return None

    def open_upstream(
        self,
        method: str,
        filename: str,
    ):
        range_header = self.headers.get("Range")
        last_error = None

        for attempt in range(
            1,
            MAX_OPEN_ATTEMPTS + 1,
        ):
            connection = http.client.HTTPSConnection(
                UPSTREAM_HOST,
                443,
                timeout=180,
                context=TLS_CONTEXT,
            )

            headers = {
                "Accept": "*/*",
                "Accept-Encoding": "identity",
                "Connection": "close",
                "User-Agent": (
                    "HG002-http1-range-proxy/1.0"
                ),
            }

            if range_header:
                headers["Range"] = range_header

            try:
                connection.request(
                    method,
                    UPSTREAM_PREFIX + filename,
                    headers=headers,
                )

                response = connection.getresponse()

                if response.status in {200, 206}:
                    return connection, response

                body = response.read(4096)

                last_error = RuntimeError(
                    f"HTTP {response.status} "
                    f"{response.reason}: "
                    f"{body[:200]!r}"
                )

                connection.close()

            except Exception as error:
                last_error = error
                connection.close()

            if attempt < MAX_OPEN_ATTEMPTS:
                delay = min(attempt * 2, 20)

                log(
                    f"Upstream open failed for "
                    f"{filename} "
                    f"({attempt}/"
                    f"{MAX_OPEN_ATTEMPTS}): "
                    f"{last_error}; "
                    f"retrying in {delay}s"
                )

                time.sleep(delay)

        raise RuntimeError(
            f"Could not open {filename} after "
            f"{MAX_OPEN_ATTEMPTS} attempts: "
            f"{last_error}"
        )

    def copy_headers(self, response) -> None:
        allowed = {
            "accept-ranges",
            "content-length",
            "content-range",
            "content-type",
            "etag",
            "last-modified",
        }

        for name, value in response.getheaders():
            if name.lower() in allowed:
                self.send_header(name, value)

    def do_HEAD(self) -> None:
        filename = self.selected_file()

        if filename is None:
            self.send_error(
                404,
                "Unknown file",
            )
            return

        connection = None

        try:
            connection, response = (
                self.open_upstream(
                    "HEAD",
                    filename,
                )
            )

            self.send_response(
                response.status,
                response.reason,
            )

            self.copy_headers(response)

            self.send_header(
                "Connection",
                "close",
            )

            self.end_headers()

        except Exception as error:
            log(
                f"HEAD failed for "
                f"{filename}: {error}"
            )

            self.send_error(
                502,
                str(error),
            )

        finally:
            self.close_connection = True

            if connection is not None:
                connection.close()

    def do_GET(self) -> None:
        filename = self.selected_file()

        if filename is None:
            self.send_error(
                404,
                "Unknown file",
            )
            return

        connection = None

        try:
            connection, response = (
                self.open_upstream(
                    "GET",
                    filename,
                )
            )

            self.send_response(
                response.status,
                response.reason,
            )

            self.copy_headers(response)

            self.send_header(
                "Connection",
                "close",
            )

            self.end_headers()

            while True:
                data = response.read(
                    BUFFER_SIZE
                )

                if not data:
                    break

                try:
                    self.wfile.write(data)
                    self.wfile.flush()

                except (
                    BrokenPipeError,
                    ConnectionResetError,
                ):
                    return

        except (
            BrokenPipeError,
            ConnectionResetError,
        ):
            return

        except Exception as error:
            log(
                f"GET failed for "
                f"{filename}: {error}"
            )

        finally:
            self.close_connection = True

            if connection is not None:
                connection.close()


class ProxyServer(ThreadingHTTPServer):
    allow_reuse_address = True
    daemon_threads = True


def main() -> None:
    server = ProxyServer(
        (LISTEN_HOST, LISTEN_PORT),
        ProxyHandler,
    )

    log(
        f"Proxy ready at "
        f"http://{LISTEN_HOST}:"
        f"{LISTEN_PORT}/"
    )

    try:
        server.serve_forever(
            poll_interval=0.5
        )

    except KeyboardInterrupt:
        pass

    finally:
        server.server_close()


if __name__ == "__main__":
    main()
