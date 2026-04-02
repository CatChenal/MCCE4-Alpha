"""
LLM interface for the agent brain (Google Gemini free tier).
"""

import re
import json
import logging
from typing import Optional
from .config import GEMINI_MODEL, get_gemini_api_key


class GeminiLLM:
    """Wrapper for Google Gemini API (free tier).

    Uses the new google-genai SDK (not the deprecated google.generativeai).
    Install: pip install google-genai
    """

    def __init__(self, api_key: Optional[str] = None):
        self.available = False
        self.client = None

        api_key = api_key or get_gemini_api_key()
        if not api_key:
            logging.info("🧠 LLM: Disabled (no GEMINI_API_KEY)")
            return

        try:
            from google import genai
            logging.getLogger("google.genai").setLevel(logging.WARNING)
            self.client = genai.Client(api_key=api_key)
            self.available = True
            logging.info(f"🧠 LLM: {GEMINI_MODEL} (free tier)")
        except ImportError:
            logging.warning("  google-genai not installed — run: pip install google-genai")
        except Exception as e:
            logging.warning(f"  Gemini init failed: {e}")

    def ask(self, prompt: str) -> Optional[str]:
        """Send prompt, return response text."""
        if not self.available:
            return None
        try:
            response = self.client.models.generate_content(
                model=GEMINI_MODEL,
                contents=prompt,
            )
            return response.text
        except Exception as e:
            logging.warning(f"  Gemini API error: {e}")
            return None

    def ask_json(self, prompt: str) -> Optional[dict]:
        """Send prompt expecting JSON response. Parses and returns dict."""
        response = self.ask(prompt)
        if not response:
            return None
        try:
            clean = response.strip()
            clean = re.sub(r'^```json\s*', '', clean)
            clean = re.sub(r'\s*```$', '', clean)
            return json.loads(clean)
        except (json.JSONDecodeError, ValueError) as e:
            logging.warning(f"  Could not parse LLM JSON: {e}")
            logging.debug(f"  Raw: {response[:300]}")
            return None
