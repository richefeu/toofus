// Copyright (C) <vincent.richefeu@3sr-grenoble.fr>
//
// This file is part of TOOFUS (TOols OFten USued)
//
// It can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note
// Without a license, the code is copyrighted by default.
// People can read the code, but they have no legal right to use it.
// To use the code, you must contact the author directly and ask permission.

#ifndef TERMINTERF_HPP
#define TERMINTERF_HPP

#include <iostream>
#include <unordered_map>
#include <variant>
#include <vector>
#include <functional>
#include <string>
#include <stack>

#include "termMenu.hpp"
#include "termPrompt.hpp"

class TermInterf {
private:
    struct Node {
        enum class Type { MENU, ACTION } type;
        std::variant<TermMenu, std::function<std::string()>> content;
    };

    struct Transition {
        std::string from;
        int choiceIndex;
        std::string next;
    };

    std::unordered_map<std::string, Node> nodes;
    std::vector<Transition> transitions;
    std::string startNode;

public:
    TermInterf() = default;

    /// Add a TermMenu node
    void addMenu(const std::string& name, const TermMenu& menu) {
        nodes[name] = Node{Node::Type::MENU, menu};
    }

    /// Add an action node (function returns name of next node or empty to return)
    void addAction(const std::string& name, const std::function<std::string()>& fn) {
        nodes[name] = Node{Node::Type::ACTION, fn};
    }

    /// Define a transition between menu items
    void link(const std::string& from, int choiceIndex, const std::string& next) {
        transitions.push_back({from, choiceIndex, next});
    }

    /// Define the starting node
    void setStart(const std::string& start) { startNode = start; }

    /// Run the terminal interface
    void run() {
        if (startNode.empty()) {
            std::cerr << "@TermInterf Error: No start node set.\n";
            return;
        }

        std::stack<std::string> history;
        std::string current = startNode;
        bool running = true;

        while (running) {
            auto it = nodes.find(current);
            if (it == nodes.end()) {
                std::cerr << "@TermInterf Error: Node '" << current << "' not found.\n";
                break;
            }

            Node& node = it->second;

            if (node.type == Node::Type::MENU) {
                TermMenu menu = std::get<TermMenu>(node.content);
                int choice = menu.run();

                bool found = false;
                for (auto& t : transitions) {
                    if (t.from == current && t.choiceIndex == choice) {
                        history.push(current);  // remember where we came from
                        current = t.next;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    std::cout << "\n(No transition for this selection, returning...)\n";
                    if (!history.empty()) {
                        current = history.top();
                        history.pop();
                    } else {
                        running = false;
                    }
                }

            } else if (node.type == Node::Type::ACTION) {
                auto& fn = std::get<std::function<std::string()>>(node.content);
                std::string next = fn();

                if (!next.empty()) {
                    current = next;
                } else if (!history.empty()) {
                    // if action returns "", go back to previous menu
                    current = history.top();
                    history.pop();
                } else {
                    running = false;
                }
            }
        }
    }
};

#endif




