#!/bin/bash

# Autor: Jorge Alberto Castro Rodríguez
# Script para utilizar bot de telegram

####==================================####
####          Bot  Telegram           ####
####==================================####

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/.env" # Cargar variables de entorno para bot de Telegram.

tg_send() {
  local msg="$1"
  curl -s -X POST "https://api.telegram.org/bot${TOKEN}/sendMessage" \
  -d "chat_id=${CHAT_ID}" \
  -d "text=${msg}" \
  -d "disable_web_page_preview=true" > /dev/null
}
