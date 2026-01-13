import sys
import logging
from pathlib import Path

# Fix python path to include src
sys.path.insert(0, str(Path.cwd() / "src"))

# Configure logging
logging.basicConfig(level=logging.INFO)

from getRPF.core.processors.rpf_extractor import ArchitectureDatabase

def test_loading():
    print("Testing architecture loading...")
    db = ArchitectureDatabase()
    
    print(f"Loaded {len(db.architectures)} architectures.")
    for arch in db.architectures:
        print(f" - {arch.protocol_name}: {len(arch.adapter_sequences)} adapters")
        
    # Verify arabidopsis is present
    arabidopsis = next((a for a in db.architectures if "arabidopsis" in a.protocol_name), None)
    if arabidopsis:
        print(f"✅ Found Arabidopsis protocol: {arabidopsis.protocol_name}")
        print(f"   Adapter: {arabidopsis.adapter_sequences[0]}")
    else:
        print("❌ Arabidopsis protocol NOT found!")

    # Verify comprehensive is present
    comp = next((a for a in db.architectures if "comprehensive" in a.protocol_name), None)
    if comp:
        print(f"✅ Found Comprehensive check with {len(comp.adapter_sequences)} adapters")
    else:
        print("❌ Comprehensive check NOT found!")

if __name__ == "__main__":
    test_loading()
