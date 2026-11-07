# Bundled Ketcher standalone editor

SEA Open embeds the official Ketcher **3.18.0** standalone application in a
same-origin iframe. The release includes the local Indigo chemistry engine.
No CDN, separate Indigo server, Node runtime, or frontend build step is needed.

- Project: https://github.com/epam/ketcher
- Release: https://github.com/epam/ketcher/releases/tag/v3.18.0
- Archive: https://github.com/epam/ketcher/releases/download/v3.18.0/ketcher-standalone-3.18.0.zip
- Archive SHA256: `484e7f10a0e74808ae5f43f0e6448db6f4d884bdd2189aa271108474167a3e7e`
- License: [Apache-2.0](LICENSE); preserve the upstream [NOTICE](NOTICE).

The archive's files are unmodified. Alternate `duo`, `popup`, and `closable`
entry points and their duplicate bundles are omitted; SEA uses `index.html`.
Retain the bundled `*.LICENSE.txt` dependency notices when redistributing.
`LICENSE` and `NOTICE` come from the matching upstream `v3.18.0` tag. These
licenses are separate from SEA Open's MIT license.

The parent interface enables only the small-molecule editor and waits for
Ketcher's `init` event before calling its asynchronous `getMolfile()` API.
RDKit validates the resulting MOL text using the same import API as files.
