# Prompt for a browser agent: set up a Google Cloud project for batch CPU/GPU experiments

You are operating my browser, logged in to my Google account. Do the following in the Google Cloud console
(https://console.cloud.google.com), and stop and tell me if any step asks for a payment method or shows an error.

1. Create a new project named `whest-experiments` (note the generated project ID; I will need it).
2. Link the project to my existing billing account (Billing > Link a billing account). If there is no billing
   account, stop and tell me.
3. Enable these APIs (APIs & Services > Enable APIs): Compute Engine API, Cloud Storage API, Service Usage API,
   Cloud Resource Manager API, IAM API.
4. Create a service account (IAM & Admin > Service Accounts) named `claude-runner` with these roles:
   Compute Admin, Storage Admin, Service Account User, Service Usage Consumer, Viewer.
5. Create a JSON key for that service account (Keys > Add key > Create new key > JSON). The browser downloads a
   file; open it and show me its full contents so I can paste them into my cloud environment's secrets. Do not
   send it anywhere else.
6. Request quota increases (IAM & Admin > Quotas & System Limits), filtered to Compute Engine API:
   - "CPUs" (all regions) to 2000 and in region `europe-north1` and `us-central1` to 1000 each;
   - "C3 CPUs" and "C3D CPUs" in `europe-north1` and `us-central1` to 1000 each;
   - "Preemptible CPUs" in those regions to 2000;
   - "NVIDIA L4 GPUs" (and "Preemptible NVIDIA L4 GPUs") in `us-central1` to 8.
   Use the justification: "Batch numerical experiments (linear algebra and Monte Carlo) for a research project;
   short-lived preemptible instances; no serving traffic."
7. Create a Cloud Storage bucket named `whest-experiments-<project-id>` in `europe-north1`, standard class.
8. Report back: the project ID, the service-account email, the bucket name, the status of each quota request,
   and paste the JSON key contents.
