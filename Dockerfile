# --- Етап 1: Збірка додатка ---
FROM gcc:latest AS builder
WORKDIR /usr/src/app
COPY . .

RUN g++ -O3 -std=c++17 NTA_lab_1.cpp -o factorizer

# --- Етап 2: Фінальний мінімальний образ ---
FROM debian:stable-slim
WORKDIR /app
COPY --from=builder /usr/src/app/factorizer .
ENTRYPOINT ["./factorizer"]